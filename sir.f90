module SIR
    use set_up_net
    ! CRITERI NODE_STATE: 0 -> SUSCEPTIBLE, 1 -> INFECTAT, 2 -> RECUPERAT
    ! CRITER REAC_RPOBS: 1-> INFECCIO, 2-> RECUPERACIO
    ! neighbors_SIR guardara la posicio del link (que especificaria 'neighbors' ) en el vector active_links
    integer, parameter :: SUSC = 0, INF = 1, RECOV = 2
    integer :: Nactive_links, Ninfected
    integer, dimension(:), allocatable :: node_state, infected_nodes, neighbors_SIR
    integer, dimension(:,:), allocatable :: active_links
    integer, parameter :: num_reacs = 2
    real(8), dimension(num_reacs) :: reac_probs, reac_rates
    real(8), dimension(0:2) :: population, pop_fraction
    real(8), dimension(0:2) :: init_pop_frac_input ! can be used or not, in the main SIR_evo call
    integer, dimension(0:2) :: init_pop_input
    
    contains
    
    subroutine read_input_SIR(file_unit)
        implicit none
        integer, intent(in) :: file_unit
        integer :: errstat, i, check_sum, still_to_add
        real(8) :: delta, lambda, init_pop_frac(0:2)
        namelist /input_SIR/ delta, lambda, init_pop_frac_input
        
        read(unit=file_unit, nml=input_SIR, iostat=errstat)
        if (errstat > 0) then
            print*, "ERROR reading namelist input_SIR from input file (code", errstat, ")"
            stop
        end if
        reac_rates = [lambda, delta]
        ! Transform the initial population fraction into number of nodes given the number of nodes in the network
        ! AIXO S'HA DE ARREGLAR:
        do i=0,2
            init_pop_input(i) = nint(dble(Nnodes) * init_pop_frac_input(i))
         !   print*, init_pop_input(i)
        enddo
        check_sum = sum(init_pop_input)
        if(check_sum.lt.Nnodes) then
            still_to_add = Nnodes - check_sum
            do while(still_to_add.gt.0)
                do i=0,2
                    if(init_pop_frac_input(i).gt.(1d-10)) then
                        init_pop_input(i) = init_pop_input(i) + 1
                        still_to_add = still_to_add - 1
                    endif
                    if(still_to_add.eq.0) exit
                enddo
            enddo
            write(*,*) "The input population fractions for SIR are not feasible. Working instead with: "
            write(*,*) "S ", dble(init_pop_input(0))/dble(Nnodes)
            write(*,*) "I ", dble(init_pop_input(1))/dble(Nnodes)
            write(*,*) "R ", dble(init_pop_input(2))/dble(Nnodes)
        endif
                    
    end subroutine read_input_SIR
    
    subroutine random_initial_pop(initial_pop,node_state,list_of_infected)
    ! Assign S, I, R population randomly according to initial pop
    ! Once each node is identified as S, I or R, active links are identified
        implicit none
        integer, intent(in) :: initial_pop(0:2)
        integer, intent(out) :: node_state(min_node:max_node),list_of_infected(min_node:max_node)
        integer :: i,j,k,l,pop_to_assign,unassigned,node_id,pos,infected
        integer, dimension(:), allocatable :: to_assign
        logical :: activate_link
        allocate(to_assign(min_node:max_node))
        unassigned = max_node
        infected = 0
        list_of_infected = 0
        do i=min_node,max_node
            to_assign(i) = i
        enddo
        do i=0,2
            pop_fraction(i) = dble(initial_pop(i))/dble(Nnodes)
            pop_to_assign = initial_pop(i)
            if(i.eq.1) Ninfected = initial_pop(i)
            do j=1,pop_to_assign
                pos = int(rand()*unassigned) + min_node
                node_id = to_assign(pos)
                node_state(node_id) = i
                if(i.eq.1) then
                    infected = infected + 1
                    list_of_infected(infected) = node_id
                endif
                ! erase bot from to_assign list
                to_assign(pos) = to_assign(unassigned)
                to_assign(unassigned) = 0
                unassigned = unassigned - 1
            enddo
        enddo
        deallocate(to_assign)
        ! Search for active links
        Nactive_links = 0
        active_links = 0
        neighbors_SIR = 0
        k = 1
        node_id = list_of_infected(k)
        do while(node_id.ne.0)
            do i=p_ini(node_id),p_fin(node_id)
                activate_link = .true.
                ! Don't activate link between two infected nodes: 
                do j=1,infected
                    if(neighbors(i).eq.list_of_infected(j)) activate_link = .false.
                enddo
                ! Don't activate link between infected and recovered:
                if(node_state(neighbors(i)).eq.2) activate_link = .false.
                ! If the neighbor is susceptible:
                if(activate_link.eqv..true.) then
                !    print*, "adding link between ", node_id, neighbors(i)
                    Nactive_links = Nactive_links + 1
                    active_links(1,Nactive_links) = node_id ! infected
                    active_links(2,Nactive_links) = neighbors(i) ! susceptible
                    neighbors_SIR(i) = Nactive_links
                    do l=p_ini(neighbors(i)),p_fin(neighbors(i))
                        if(neighbors(l).eq.node_id) neighbors_SIR(l) = Nactive_links
                    enddo
                endif
            enddo
            k = k + 1
            node_id = list_of_infected(k)
        enddo
    end subroutine random_initial_pop
    
    subroutine compute_reac_probs()
        ! Es faria servir en cas que la probabilitat de reaccio depengui del numero de nodes infectats.
        ! Si les probabilitats son directament les rates i son constants, no caldra cridar aquesta subrutina
        ! No es normalizen les probs (no cal fer-ho pel tower sampling de triar la reaccio) // pero si cal per trobar tau!
        implicit none
        real(8) :: sum_aux
        integer :: i
        reac_probs(1) = dble(Nactive_links) * reac_rates(1)
        reac_probs(2) = dble(Ninfected) * reac_rates(2)
        sum_aux = sum(reac_probs)
        reac_probs = reac_probs/sum_aux
    end subroutine compute_reac_probs
    
    real(8) function compute_reac_time(a0)
        !use mtmod
        implicit none
        real(8), intent(in) :: a0
        !compute_reac_time = - log(grnd())/a0
        compute_reac_time = - log(rand())/a0
    end function compute_reac_time
    
    integer function choose_reac(sum_probs,reac_probs)
    ! Returns the integer identifying the reaction: INFECTION (1) or RECOVERING (2)
        !use mtmod
        implicit none
        integer i
        real(8), intent(in) :: sum_probs
        real(8), dimension(num_reacs), intent(in) :: reac_probs
        real(8) alea,sum_aux
        ! Recerca lineal:
        sum_aux = 0d0
        !alea = grnd()*sum_probs
        alea = (1d0-rand())*sum_probs
        !print*, "alea", alea
        search: do i=1,num_reacs
            sum_aux = sum_aux + reac_probs(i)
            if(alea.le.sum_aux) then
                choose_reac = i
                exit search
            endif
        enddo search
    end function choose_reac

    subroutine SIR_evolution(initial_pop,file_unit,max_inf,max_infrec)
        implicit none
        integer, intent(in) :: initial_pop(0:2), file_unit
        real(8), intent(out) :: max_inf, max_infrec
        ! internal:
        integer :: i,j,reac_i,max_iters
        real(8) :: sum_probs, tau, time
        
        allocate(node_state(min_node:max_node))
        allocate(infected_nodes(min_node:max_node))
        allocate(active_links(2,num_links))
        allocate(neighbors_SIR(2*num_links))
        
        population = initial_pop
        call random_initial_pop(initial_pop,node_state,infected_nodes)
        
        max_inf = 0d0
        max_infrec = 0d0
        !print*, reac_rates
        !print*, "node_state", node_state
        !print*, "infected_nodes", infected_nodes
        !print*, " "
        !print*, "Active links (inf)", active_links(1,:)
        !print*, "Active links (susc)", active_links(2,:)
        !print*
        !print*, "n", neighbors
        !print*, "n_SIR", neighbors_SIR
        max_iters = Nnodes * 1000
        time = 0
        do i=1,max_iters
            call compute_reac_probs()
            !print*, "reac_probs", reac_probs(:), reac_rates(1)
            sum_probs = sum(reac_probs)
            !print*, "sum_probs", sum_probs
            tau = compute_reac_time(sum_probs)
            !print*, "tau", tau
            reac_i = choose_reac(sum_probs,reac_probs)
         !   print*, "chosen reaction", reac_i
            call update_system(reac_i)
            !print*, " "
         !   print*, "UPDATE AT ITER", i, " ******************************************************** "
         !   print*, "node_state", node_state
         !   print*, "infected_nodes", infected_nodes
         !   print*, " "
         !   print*, "Active links (inf)", active_links(1,:)
         !   print*, "Active links (susc)", active_links(2,:)
         !   print*
         !   print*, "n", neighbors
         !   print*, "n_SIR", neighbors_SIR
            time = time + tau
            ! Chex maximal populations
            if(pop_fraction(1).gt.max_inf) max_inf = pop_fraction(1)
            if(sum(pop_fraction(1:2)).gt.max_infrec) max_infrec = sum(pop_fraction(1:2))
            if(file_unit.ne.1000) then
                write(file_unit,*) i, time, pop_fraction(:)
            endif
            if(population(RECOV).eq.Nnodes) then
                !write(*,*) " ****** All population recovered, exiting SIR at iter ", i, " ****** "
                exit
            endif
            if(population(INF).eq.0) then
                !write(*,*) " ****** 0 population infected, exiting SIR at iter ", i, " ****** "
                exit
            endif
        enddo
        if(allocated(node_state)) deallocate(node_state)
        if(allocated(infected_nodes)) deallocate(infected_nodes)
        if(allocated(active_links)) deallocate(active_links)
        if(allocated(neighbors_SIR)) deallocate(neighbors_SIR)
    end subroutine SIR_evolution
    
    subroutine update_system(reac_i)
    ! actualitzar pop_fractions !!
        implicit none
        integer, intent(in) :: reac_i
        integer :: i,j
        integer :: infection_link, infector_node, infected_node
        integer :: pos_rec_node, recovered_node, susc_node
        
        select case(reac_i)
            case(1) ! INFECTION
                ! Triar un node infectat, despres infectar un dels seus veins susceptibles -> Triar un link actiu
                ! Un cop fet, actualitzar llista infectats (sumar node infectat) i llista links actius (treure link actiu i afegir nous)
                infection_link = int(rand()*Nactive_links) + 1 
                infector_node = active_links(1,infection_link)
                infected_node = active_links(2,infection_link)
                !print*, "Node ", infector_node, " is infecting node ", infected_node
                if(infected_node.eq.0) then
                    !print*, "ac_links", active_links
                    print*
                    print*
                    print*, "inf_nodes", infected_nodes
                    print*, "Nactive_links", Nactive_links
!                    print*, reac_probs
                endif
                node_state(infected_node) = INF
                ! Add infected node to the list of infected nodes:
                Ninfected = Ninfected + 1
                infected_nodes(Ninfected) = infected_node
                do i=p_ini(infected_node),p_fin(infected_node)
                    if(neighbors_SIR(i).ne.0) then
                        ! Update position of swaped link (the last one) for the position to wich is going (the erased link position)
                        if(neighbors_SIR(i).ne.Nactive_links) then
                            do j=p_ini(active_links(1,Nactive_links)),p_fin(active_links(1,Nactive_links))
                                if(neighbors(j).eq.active_links(2,Nactive_links)) neighbors_SIR(j) = neighbors_SIR(i)
                            enddo
                            do j=p_ini(active_links(2,Nactive_links)),p_fin(active_links(2,Nactive_links))
                                if(neighbors(j).eq.active_links(1,Nactive_links)) neighbors_SIR(j) = neighbors_SIR(i)
                            enddo
                        endif
                        do j=p_ini(active_links(1,neighbors_SIR(i))),p_fin(active_links(1,neighbors_SIR(i)))
                            if(neighbors(j).eq.infected_node) then
                                !print*, " ----------- erasing connect. bt. ",active_link(1,neighbors_SIR(i)), " and ", infected_node, "in position ", j , " of n_SIR"
                                neighbors_SIR(j) = 0
                            endif
                        enddo
                        active_links(:,neighbors_SIR(i)) = active_links(:,Nactive_links)
                        active_links(:,Nactive_links) = 0
                        Nactive_links = Nactive_links - 1
                        neighbors_SIR(i) = 0
                    endif
                enddo
                ! Erase active link going from infector to infected form the links position list:
                do i=p_ini(infector_node),p_fin(infector_node)
                    if(neighbors(i).eq.infected_node) neighbors_SIR(i) = 0
                enddo
                ! Search new links:
                do i=p_ini(infected_node),p_fin(infected_node)
                    if(node_state(neighbors(i)).eq.SUSC) then
                        Nactive_links = Nactive_links + 1
                        active_links(1,Nactive_links) = infected_node
                        active_links(2,Nactive_links) = neighbors(i)
                        neighbors_SIR(i) = Nactive_links
                        do j=p_ini(neighbors(i)),p_fin(neighbors(i))
                            if(neighbors(j).eq.infected_node) neighbors_SIR(j) = Nactive_links
                        enddo
                    endif
                enddo
                ! Update population fractions:
                pop_fraction(SUSC) = pop_fraction(SUSC) - 1d0/dble(Nnodes)
                pop_fraction(INF) = pop_fraction(INF) + 1d0/dble(Nnodes)
                population(SUSC) = population(SUSC) - 1
                population(INF) = population(INF) + 1
            case(2) ! RECOVERY
                ! Choose an infected node:
                pos_rec_node = int(rand()*Ninfected) + min_node
                recovered_node = infected_nodes(pos_rec_node)
                node_state(recovered_node) = RECOV
             !   print*, "node ", recovered_node, " is recovering "
                ! Erase recovered node from list of infected nodes:
                infected_nodes(pos_rec_node) = infected_nodes(Ninfected)
                infected_nodes(Ninfected) = 0
                Ninfected = Ninfected - 1
                ! Erase active links that stem from recovered node:
                if(Nactive_links.gt.0) then
                    do i=p_ini(recovered_node),p_fin(recovered_node)
                        !print*, "neighbor: ", neighbors(i)
                        !print*, "position in active links: ", neighbors_SIR(i)
                        if(neighbors_SIR(i).ne.0) then
                            ! Update position of swaped link (the last one) for the position to which is going (the erased link position)
                            !print*, "active link infected ", active_links(1,neighbors_SIR(i))
                            !print*, "active link susceptible ", active_links(2,neighbors_SIR(i))
                            susc_node = active_links(2,neighbors_SIR(i))
                            if(neighbors_SIR(i).ne.Nactive_links) then
                                do j=p_ini(active_links(1,Nactive_links)),p_fin(active_links(1,Nactive_links))
                                    if(neighbors(j).eq.active_links(2,Nactive_links)) neighbors_SIR(j) = neighbors_SIR(i)
                                enddo
                                do j=p_ini(active_links(2,Nactive_links)),p_fin(active_links(2,Nactive_links))
                                    if(neighbors(j).eq.active_links(1,Nactive_links)) neighbors_SIR(j) = neighbors_SIR(i)
                                enddo
                            endif
                            active_links(:,neighbors_SIR(i)) = active_links(:,Nactive_links)
                            active_links(:,Nactive_links) = 0
                            Nactive_links = Nactive_links - 1
                            do j=p_ini(susc_node),p_fin(susc_node)
                                if(neighbors(j).eq.recovered_node) neighbors_SIR(j) = 0
                            enddo
                            do j=p_ini(recovered_node),p_fin(recovered_node)
                                if(neighbors(j).eq.susc_node) neighbors_SIR(j) = 0
                            enddo
                        endif
                    enddo
                endif
                ! Update population fractions:
                pop_fraction(INF) = pop_fraction(INF) - 1d0/dble(Nnodes)
                pop_fraction(RECOV) = pop_fraction(RECOV) + 1d0/dble(Nnodes)
                population(INF) = population(INF) - 1
                population(RECOV) = population(RECOV) + 1
            end select
        end subroutine update_system
                
        
end module SIR
        
