program main
    use set_up_net
    use net_props
    use SIR
    implicit none
    character(len=50) :: input_name
    real(8) :: cbar
    integer :: initial_pop_SIR(0:2)
    
    ! To execute the program >> main.x input_file. Otherwise, an error will occur.
    if (command_argument_count() == 0) stop "ERROR: Cridar fent >> ./main.x input_net.dat"
    call get_command_argument(1, input_name)
      
    open(10, file=input_name)
    call read_net(10)
    close(10)
    
   ! call comp_degree_distros()
   ! write(*,*)
   ! write(*,*) " Degree distr: ", degree_distr
   ! write(*,*) " Cumulative degree distr: ", c_degree_distr
   ! write(*,*) " Complementary cumulative degree distr: ", cc_degree_distr
    
   ! call comp_avg_nn_degree_distro()
   ! write(*,*)
   ! write(*,*) " Average nearest neighbor degree distro: ", knn
    
   ! call count_triangles()
    
   ! call comp_clustering()
   ! cbar = avg_clustering()
   ! print*, "cbar ", cbar
    
    
    ! SIR evolution: 
    
    ! aixo s'hauria de posar en un input namelist o algo
    initial_pop_SIR = [2, 3, 0]
    call srand(1234)
    
    open(10, file="input_SIR.txt")
    call read_input_SIR(10)
    close(10)
    call SIR_evolution(initial_pop_SIR)
    
end program
