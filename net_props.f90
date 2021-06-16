module net_props
    use set_up_net
    implicit none
    integer, dimension(:), allocatable :: N_nodes_degree
    real(8), dimension(:), allocatable :: degree_distr, c_degree_distr, cc_degree_distr, knn, clustering
    
    contains
    
    subroutine comp_degree_distros()
        ! Computes the degree distribution, cumulative degree distribution, complementary cumulative degree distribution
        implicit none
        integer :: i, count_degree
        
        allocate(N_nodes_degree(min_degree:max_degree))
        allocate(degree_distr(min_degree:max_degree))
        allocate(c_degree_distr(min_degree:max_degree))
        allocate(cc_degree_distr(min_degree:max_degree))

        degree_distr = 0d0
        do i=min_node,max_node
            degree_distr(degree(i)) = degree_distr(degree(i)) + 1d0
            N_nodes_degree(degree(i)) = N_nodes_degree(degree(i)) + 1
        enddo
        degree_distr = degree_distr/sum(degree_distr)
        
        do i=min_degree,max_degree
            if(i.eq.min_degree) then
                c_degree_distr(i) = degree_distr(i)
                cc_degree_distr(i) = 1d0
            else
                c_degree_distr(i) = c_degree_distr(i-1) + degree_distr(i)
                cc_degree_distr(i) = cc_degree_distr(i-1) - degree_distr(i)
            endif
        enddo
    end subroutine comp_degree_distros
    
    subroutine comp_avg_nn_degree_distro()
         ! Computes the average nearest neighbor degree distribution, knn
         implicit none
         integer :: i, j
         
         allocate(knn(min_degree:max_degree))
         knn = 0d0
         
         do i=min_node,max_node
             do j=p_ini(i),p_fin(i)
                 knn(degree(i)) = knn(degree(i)) + degree(neighbors(j))/(dble(N_nodes_degree(degree(i))*degree(i)))
             enddo
         enddo
     end subroutine comp_avg_nn_degree_distro
     
     subroutine count_triangles()
         ! te en compte que neighbors esta ordenat, i.e. neigbors de 1 -> 2 3 4, NO 1 -> 2 4 3 !!!
         implicit none
         integer :: i, j, k, l, triangles, apex_1, apex_2
         
         triangles = 0
         do i=min_node,max_node
         !    do j=p_ini(i),p_fin(i)-1 ! no se benbe pq ho faig, pero funciona per contar triangles. tmb funciona si no ho faig
              do j=p_ini(i),p_fin(i)
                 apex_1 = neighbors(j)
                 if(apex_1.gt.i) then 
                     do k=p_ini(apex_1),p_fin(apex_1)
                         !do l=j,p_fin(i)
                         apex_2 = neighbors(k)
                         if(apex_2.gt.apex_1) then ! ckeck they have i as common neighbor
                             do l=p_ini(apex_2),p_fin(apex_2)
                                 if(neighbors(l).eq.i) then
                                     triangles = triangles + 1
                                    ! print*, "triangle by", i, apex_1, apex_2
                                 endif
                             enddo
                         endif
                     enddo
                 endif
             enddo
         enddo
         print*, "triangles ", triangles
     end subroutine count_triangles
     
     subroutine comp_clustering()
     ! te en compte que neighbors esta ordenat, i.e. neigbors de 1 -> 2 3 4, NO 1 -> 2 4 3 !!!
         implicit none
         integer :: i, j, k, l, triangles_node_i, apex_1, apex_2
         
         allocate(clustering(min_degree:max_degree))
         clustering = 0d0
         
         do i=min_node,max_node
             triangles_node_i = 0
              do j=p_ini(i),p_fin(i)
                 apex_1 = neighbors(j)
        !         if(apex_1.gt.i) then 
                     do k=p_ini(apex_1),p_fin(apex_1)
                         !do l=j,p_fin(i)
                         apex_2 = neighbors(k)
        !                 if(apex_2.gt.apex_1) then ! ckeck they have i as common neighbor
                             do l=p_ini(apex_2),p_fin(apex_2)
                                 if(neighbors(l).eq.i) then
                                     triangles_node_i = triangles_node_i + 1
                                 endif
                             enddo
        !                 endif
                     enddo
        !         endif
             enddo
             !if(mod(triangles_node_i,2).ne.0) print*, "aiba!", triangles_node_i
             triangles_node_i = triangles_node_i/2
             !print*, "node ", i, "triangles ", triangles_node_i
             if(degree(i).gt.1) then
                 clustering(degree(i)) = clustering(degree(i)) + 2d0*dble(triangles_node_i)/&
                 (dble(degree(i)*(degree(i)-1))*degree_distr(degree(i))*dble(Nnodes))
             endif
         enddo
    end subroutine comp_clustering
    
    real(8) function avg_clustering()
        implicit none
        integer :: i
        real(8) :: alt,k_sq_avg, k_avg
        avg_clustering = 0d0
        do i=min_degree,max_degree
            avg_clustering = avg_clustering + degree_distr(i)*clustering(i)
        enddo
        k_sq_avg = 0d0
        k_avg = 0d0
        do i=min_node,max_node
            k_sq_avg = k_sq_avg + degree(i)**2/dble(Nnodes)
            k_avg = k_avg + degree(i)/dble(Nnodes)
        enddo
        alt = (k_sq_avg - k_avg)**2/k_avg**3
        alt = alt/dble(Nnodes)
        print*, "altrenate cbar ", alt
        return
   end function avg_clustering
     
end module net_props
