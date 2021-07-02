program main
    use set_up_net
    use net_props
    use SIR
    implicit none
    character(len=50) :: input_name, seed_char
    real(8) :: cbar, avg_max_inf, avg_max_infrec, max_inf, max_infrec, lambda_max, lambda_min, dlambda, dlambda_lt1, lambda
    integer :: i,j, initial_pop_SIR(0:2), seed, Nrea, rea_write, iter_lambda, iter_lambda_lt1, lambda_int, file_unit
    character(len=4) :: lambda_tag
    
    ! To execute the program >> main.x input_file. Otherwise, an error will occur.
    if (command_argument_count().ne.2) stop "ERROR: Cridar fent >> ./main.x input_net.dat SEED"
    call get_command_argument(1, input_name)
    call get_command_argument(2, seed_char)
    read(seed_char, *) seed
      
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
    !initial_pop_SIR = [4, 1, 0]
!    call srand(8317909)
!    open(10, file="input_SIR.txt")
!    call read_input_SIR(10)
!    close(10)
!    reac_rates = [3.6d0, 1d0]
!    open(10, file="pop_frac_evo.dat")
!    call SIR_evolution(init_pop_input,10,max_inf,max_infrec)
!    close(10)
    
    ! Evolució temporal, differents rates
    lambda_min = 0.1d0
    lambda_max = 5d0
    dlambda = 0.25d0
    dlambda_lt1 = 0.1d0
    iter_lambda_lt1 = int((1d0-0.1d0)/dlambda_lt1)
    iter_lambda = int((lambda_max - 1d0)/dlambda)
    ! Nrea = 100 ! if computing avg max infected starting from an infected node
    Nrea = 5
    rea_write = 5
  
    open(10, file="input_SIR.txt")
    call read_input_SIR(10)
    close(10)
    
    open(11, file="max_inf_infrec_lambda.dat")
    lambda = 0d0
    do j=0,iter_lambda_lt1
        lambda = lambda_min + j*dlambda_lt1
        ! Overwrite reaction rate from input_file:
        reac_rates(1) = lambda
        avg_max_inf = 0
        avg_max_infrec = 0
        !open(10, ...) ! ??? amb el nom de lambda
        lambda_int = lambda*10
        if(mod(lambda_int,1).eq.0) then
            write(lambda_tag,'(I0)') lambda_int
            file_unit = 10
            open(file_unit, file="pop_frac_evo_lambda_"//trim(lambda_tag)//".dat")
        else
            file_unit = 1000
        endif
        do i=1,Nrea
            print*, "seed", seed, "lambda", lambda
            call srand(seed)
            if(i.le.rea_write) then ! make simulation write temporal evolution output
                call SIR_evolution(init_pop_input,file_unit,max_inf,max_infrec)
                if(mod(lambda_int,1).eq.0) then
                    write(file_unit,*)
                    write(file_unit,*)
                endif
            else
                call SIR_evolution(init_pop_input,1000,max_inf,max_infrec)
            endif
            avg_max_inf = avg_max_inf + max_inf
            avg_max_infrec = avg_max_infrec + max_infrec
            seed = int(10000000*rand())
        enddo
        if(mod(lambda_int,5).eq.0) then
            close(file_unit)
        endif
        avg_max_inf = avg_max_inf/dble(Nrea)
        avg_max_infrec = avg_max_infrec/dble(Nrea)
        write(11,*) lambda, avg_max_inf, avg_max_infrec
    enddo
    
    do j=1,iter_lambda
        lambda = 1d0 + j*dlambda
        ! Overwrite reaction rate from input_file:
        reac_rates(1) = lambda
        avg_max_inf = 0
        avg_max_infrec = 0
        lambda_int = lambda*10
        if(mod(lambda_int,5).eq.0) then
            write(lambda_tag,'(I0)') lambda_int
            file_unit = 10
            open(file_unit, file="pop_frac_evo_lambda_"//trim(lambda_tag)//".dat")
        else
            file_unit = 1000
        endif
        do i=1,Nrea
            print*, "seed", seed, "lambda", lambda
            call srand(seed)
            if(i.le.rea_write) then ! make simulation write temporal evolution output
                call SIR_evolution(init_pop_input,file_unit,max_inf,max_infrec)
                if(mod(lambda_int,5).eq.0) then
                    write(file_unit,*)
                    write(file_unit,*)
                endif
            else
                call SIR_evolution(init_pop_input,1000,max_inf,max_infrec)
            endif
            avg_max_inf = avg_max_inf + max_inf
            avg_max_infrec = avg_max_infrec + max_infrec
            seed = int(10000000*rand())
        enddo
        if(mod(lambda_int,5).eq.0) then
            close(file_unit)
        endif
        avg_max_inf = avg_max_inf/dble(Nrea)
        avg_max_infrec = avg_max_infrec/dble(Nrea)
        write(11,*) lambda, avg_max_inf, avg_max_infrec
    enddo
    
    
end program
