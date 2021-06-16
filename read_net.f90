module set_up_net
    implicit none
    integer :: num_links, min_node, max_node, min_degree, max_degree, Nnodes
    integer, dimension(:), allocatable :: p_ini, p_fin, degree, neighbors, nodes
    
    contains
    
    subroutine read_net(nfi)
        ! reads the net in the input file and sets up global variables, defined in this module
        implicit none
        integer, intent(in) :: nfi ! net file index
        integer :: read_stat,i,j,count_links,count_links2,aux_int,aux_int2,swaps
        integer, dimension(:,:), allocatable :: check_reps
        integer, dimension(2) :: first_link
        logical :: reiterate_links
        
        ! Assume that links are not reiterated, if found one, then work as reiterated:
        reiterate_links = .false.
        first_link = [0, 0]
        
        ! 1st lecture: define max_node and number of links (count_links)
        max_node = 0
        min_node = 10
        count_links = 0
        do while(1.EQ.1)
            read(nfi,*,iostat=read_stat) i,j
            if((first_link(1).eq.0).and.(first_link(2).eq.0).and.(reiterate_links.eqv..false.)) then !will stop checking when rei_links=t
               first_link = [i,j]
            else
                if((first_link(1).eq.j).and.(first_link(2).eq.i)) reiterate_links = .true.
            endif
            if(IS_IOSTAT_END(read_stat)) then
                exit
            else
                count_links = count_links + 1
                aux_int = max(i,j)
                aux_int2 = min(i,j)
                if(aux_int.gt.max_node) max_node = aux_int
                if(aux_int2.lt.min_node) min_node = aux_int2
            endif
        enddo
    
        ! ASSUME NO REPETITIONS, THIS IS VERY DEMANING FOR LARGE NETWORKS
        ! Now that the number of links is known, check if there is any repeated link in net, and write an auxiliar net_no_repeat.dat 
        ! that will be used during the rest of the execution.
!        rewind(nfi)
!        allocate(check_reps(count_links,2))
!        count_links2 = 0
        ! 2nd lecture: check for repetitions
!        do while(1.EQ.1)
!            read(nfi,*,iostat=read_stat) i,j
!            if(IS_IOSTAT_END(read_stat)) then
!                exit
!            else
!                count_links2 = count_links2 + 1
!                check_reps(count_links2,:) = [i, j]
!                if(count_links2.gt.1) then ! check for repetitions
!                    aux_int = count_links2-1
!                    loop_reps: do i=1,aux_int
!                        if((check_reps(count_links2,1).eq.check_reps(i,1)).and.(check_reps(count_links2,2).eq.check_reps(i,2))) then
!                            count_links2 = count_links2-1
!                            exit loop_reps
!                        endif
!                    enddo loop_reps
!                endif
!            endif
!        enddo
    
        ! Process repetitions
!        if(count_links2.lt.count_links) then
!            count_links = count_links2
!            close(nfi)
!            open(nfi, file="net_no_repeat.dat")
!            do i=1,count_links
!                write(nfi,*) check_reps(i,:)
!            enddo
!            deallocate(check_reps)
!            rewind(nfi) ! and use net_no_repeat.dat from now on
!        else
!            deallocate(check_reps)
!            rewind(nfi) ! and keep using net.dat
!        endif
    
        rewind(nfi) ! al haver commentat la part de check reps
        allocate(degree(min_node:max_node))
        allocate(nodes(min_node:max_node)) ! variable per saber si tots els numeros entre min i max node es fan servir realment
        nodes = 0
        print*, "max node: ", max_node
        print*, "min node: ", min_node
        if(reiterate_links.eqv..true.) then
            print*, "links: ", count_links/2
        else
            print*, "links: ", count_links
        endif
    
        degree = 0
        do while(1.eq.1)
            read(nfi,*,iostat=read_stat) i,j
            if(IS_IOSTAT_END(read_stat)) then
                exit
            else
                degree(i) = degree(i) + 1
                if(nodes(i).eq.0) nodes(i) = 1
                if(reiterate_links.eqv..false.) then
                    degree(j) = degree(j) + 1 !fer-ho si els links no es repeteixen, no fer-ho si a net.dat hi ha e.g. 1 2 i despres 2 1
                    if(nodes(j).eq.0) nodes(j) = 1
                endif
            endif
        enddo
    
        Nnodes = sum(nodes)
        print*, "number of nodes", Nnodes
        deallocate(nodes)
        min_degree = minval(degree)
        max_degree = maxval(degree)
        rewind(nfi)
        !print*, "degree", degree(:)
        if(reiterate_links.eqv..true.) then
            allocate(neighbors(count_links))
        else
            allocate(neighbors(2*count_links))
        endif
        neighbors = 0
        allocate(p_ini(min_node:max_node))
        allocate(p_fin(min_node:max_node))
    
        ! Assign p_ini:
        p_ini(min_node) = 1
        do i=min_node+1,max_node
            p_ini(i) = p_ini(i-1) + degree(i-1)
            p_fin(i-1) = p_ini(i-1) - 1
        enddo
        p_fin(max_node) = p_ini(max_node) - 1
        ! Assign p_fin and neighbors:
        do while(1.eq.1)
            read(nfi,*,iostat=read_stat) i,j
            if(IS_IOSTAT_END(read_stat)) then
                exit
            else
                if(reiterate_links.eqv..true.) then
                    if(j.gt.i) then
                        p_fin(i) = p_fin(i) + 1
                        p_fin(j) = p_fin(j) + 1
                        neighbors(p_fin(i)) = j
                        neighbors(p_fin(j)) = i
                    endif
                else
                    p_fin(i) = p_fin(i) + 1
                    p_fin(j) = P_fin(j) + 1
                    neighbors(p_fin(i)) = j
                    neighbors(p_fin(j)) = i
                endif
            endif
        enddo
    
       ! print*, "neighbors", neighbors(:)
       ! print*, "p_ini", p_ini(:)
       ! print*, "p_fin", p_fin(:)
        
        ! Ordena neighbors de mes petit a mes gran per cada node
        do i=min_node,max_node
            swaps = -1
            do while(swaps.ne.0)
                swaps = 0
                do j=p_ini(i),p_fin(i)-1
                    if(neighbors(j).gt.neighbors(j+1)) then
                        aux_int = neighbors(j+1)
                        neighbors(j+1) = neighbors(j)
                        neighbors(j) = aux_int
                        swaps = swaps + 1
                    endif
                enddo
            enddo
        enddo
        
       ! print*, neighbors(:)
        
    end subroutine read_net
end module set_up_net
