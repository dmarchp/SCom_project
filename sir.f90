module SIR
    use set_up_net
    ! CRITERI NODE_STATE: 0 -> SUSCEPTIBLE, 1 -> INFECTAT, 2 -> RECUPERAT
    ! CRITER REAC_RPOBS: 1-> INFECCIO, 2-> RECUPERACIO
    integer :: Nactive_links
    integer, dimension(:), allocatable :: node_state, infected_nodes
    integer, dimension(:,:), allocatable :: active_links
    integer, parameter :: num_reacs = 2
    real(8), dimension(num_reacs) :: reac_probs, reac_rates
    real(8), dimension(0:2) :: population, pop_fraction
    real(8), dimension(0:2) :: init_pop_frac_input ! can be used or not, in the main SIR_evo call
    
    contains
    
    subroutine read_input_SIR(file_unit)
        implicit none
        integer, intent(in) :: file_unit
        integer :: errstat, i
        real(8) :: delta, lambda, init_pop_frac(0:2)
        namelist /input_SIR/ delta, lambda, init_pop_frac
        
        read(unit=file_unit, nml=input_SIR, iostat=errstat)
        if (errstat > 0) then
            print*, "ERROR reading namelist input_SIR from input file (code", errstat, ")"
            stop
        end if
        reac_rates = [lambda, delta]
        ! Transform the initial population fraction into number of nodes given the number of nodes in the network
        ! AIXO S'HA DE ARREGLAR:
        do i=0,2
            init_pop_frac_input(i) = int(dble(Nnodes) * init_pop_frac_input(i))
        enddo
    end subroutine read_input_SIR
    
    subroutine random_initial_pop(initial_pop,node_state,list_of_infected)
    ! Assign S, I, R population randomly according to initial pop
    ! Once each node is identified as S, I or R, active links are identified
        implicit none
        integer, intent(in) :: initial_pop(0:2)
        integer, intent(out) :: node_state(min_node:max_node),list_of_infected(min_node:max_node)
        integer :: i,j,k,pop_to_assign,unassigned,node_id,pos,infected
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
            pop_to_assign = initial_pop(i)
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
        k = 1
        node_id = list_of_infected(k)
        do while(node_id.ne.0)
            do i=p_ini(node_id),p_fin(node_id)
                activate_link = .true.
                do j=1,infected
                    if(neighbors(i).eq.list_of_infected(j)) activate_link = .false.
                enddo
                if(activate_link.eqv..true.) then
                    print*, "adding link between ", node_id, neighbors(i)
                    Nactive_links = Nactive_links + 1
                    active_links(1,Nactive_links) = node_id ! infected
                    active_links(2,Nactive_links) = neighbors(i) ! susceptible
                endif
            enddo
            k = k + 1
            node_id = list_of_infected(k)
        enddo
    end subroutine random_initial_pop
    
    subroutine compute_reac_probs()
        ! Es faria servir en cas que la probabilitat de reaccio depengui del numero de nodes infectats.
        ! Si les probabilitats son directament les rates i son constants, no caldra cridar aquesta subrutina
        ! No es normalizen les probs (no cal fer-ho pel tower sampling de triar la reaccio)
        implicit none
        integer :: i
        do i=1,num_reacs
            reac_probs = pop_fraction(i) * reac_rates(i)
        enddo
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
        alea = rand()*sum_probs
        search: do i=1,num_reacs
            sum_aux = sum_aux + reac_probs(i)
            if(alea.le.sum_probs) then
                choose_reac = i
                exit search
            endif
        enddo search
    end function choose_reac

    subroutine SIR_evolution(initial_pop)
        implicit none
        integer, intent(in) :: initial_pop(0:2)
        ! internal:
        integer :: i,j,reac_i
        real(8) :: sum_probs, tau
        
        allocate(node_state(min_node:max_node))
        allocate(infected_nodes(min_node:max_node))
        allocate(active_links(2,num_links))
        
        population = initial_pop
        call random_initial_pop(initial_pop,node_state,infected_nodes)
        
        print*, node_state
        print*, infected_nodes
        print*, " "
        print*, active_links(1,:)
        print*, active_links(2,:)
        !do i=1,max_iters
        !    call compute_reac_probs()
        !    sum_probs = sum(reac_probs)
        !    tau = compute_reac_time(sum_probs)
        !    reac_i = choose_reac(a0,reac_probs)
        !    call update_system(reac_i,tau, ? )! TO IMPLEMENT: AQUI S'HAURA DE CRIDAR ALS SWAPS DELS VECTORS infected_nodes, active_links
        !enddo
    end subroutine SIR_evolution
        
end module SIR
        
