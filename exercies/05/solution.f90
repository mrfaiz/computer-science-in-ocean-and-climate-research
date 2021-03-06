program PreditorPreySpace
    use precision
    USE OMP_LIB
    implicit none

    integer(4) :: N, steps, calculate_steps,i,j
    real :: start_time, end_time
    real(kind=wp) :: alpha,beta,gamma,delta,lambda,mu,kappa,t0,T,dt,delta_single_box_length,delta_time,diffusion_factor
    
    real(kind=wp), allocatable :: boxes(:),times(:),D(:,:),preditor_matrix(:,:), prey_matrix(:,:)
    real(kind=wp), allocatable :: xk(:), yk(:),prey_f1(:),preditor_f2(:)

    namelist /model_parameters/ alpha, beta, gamma, delta, lambda, mu
    namelist /spatial_parameters/ kappa, N
    namelist /time_parameters/ t0, T, dt

    call cpu_time(start_time)
    open(unit = 20, file = 'predatorprey.nml', action = 'read')
    read(20,nml=model_parameters)
    read(20,nml=spatial_parameters)
    read(20,nml=time_parameters)

    steps = int((T-t0)/dt)
   
    ! Make sure in your code that the condition n ≥ 2N 2 κT is always satisfied
    calculate_steps = int(2*N*N*kappa*T)

    if(steps<calculate_steps) then
        steps = calculate_steps
    endif
    
    print*,steps
    
  ! Allocation
    allocate(boxes(N))
    allocate(times(steps))
    allocate(D(N,N))
    allocate(preditor_matrix(n,N))
    allocate(prey_matrix(steps,N))
    allocate(yk(N))
    allocate(xk(N))
    allocate(preditor_f2(N))
    allocate(prey_f1(N))
    ! Box Data
    delta_single_box_length = 1.0/DFLOAT(N) ! also h
    do i = 1, N
        boxes(i) = i * delta_single_box_length
    enddo
    
    ! Calculate Delta_t and populate data to time vector
    delta_time = (T - 0.0)/ DFLOAT(n)
    
    ! Calculate Diffusion multiplication Factor
    diffusion_factor = kappa/(delta_single_box_length * delta_single_box_length)
    ! print*, diffusion_factor
   
    ! populate data in diffussion matrix
    do i = 1, N ! row
        do j = 1 , N !colum
            if (i==j) then
                if (i == 1 .OR. i == N) then
                    D(i,j) = (-1)*diffusion_factor
                else
                    D(i,j) = (-2)*diffusion_factor
                endif
            else
                if(i == j+1 .OR. j == i+1)then
                    D(i,j) = 1 * diffusion_factor
                else
                    D(i,j) = 0.0
                endif
            endif
        enddo
    enddo

    ! populate initial(random data) data inside preditor, prey and times matrix

    do i = 1, N 
        if(i<=(N/2)) then
            preditor_matrix(1,i) = 0.2
            prey_matrix(1,i) = 1.0
        else
            preditor_matrix(1,i) = 2.0
            prey_matrix(1,i) = 0.1
        endif
    enddo

    times(1) = 0.0 

    ! Main Loop


    !  for parallalizations 

    !$OMP DO  
    do i = 2, n 
        xk = prey_matrix(i-1,:)
        yk = preditor_matrix(i-1,:)

        ! f 1 (x k , y k ) = Dx k + x k ∗ (α − βy k − λx k )
        ! f 2 (x k , y k ) = Dy k + y k ∗ (δx k − γ − μy k )
        prey_f1 = matmul(D,xk) + xk * (alpha - (beta * yk) - (lambda * xk))
        preditor_f2 = matmul(D,yk) + yk * ((delta * xk) - gamma - (mu * yk))  

        prey_matrix(i,:) = xk + delta_time * prey_f1
        preditor_matrix(i,:) = yk + delta_time * preditor_f2

        ! tk+1 = tk + delta_t * k
        times(i) = times(i-1) + delta_time
    enddo
    !$OMP END DO

    close(20)
    call cpu_time(end_time)
    print '("Time = ",f6.3," seconds.")', end_time - start_time
end program PreditorPreySpace

! For double precision 
!   
!   Without OMP
    !   Time =  3.147 seconds.
!   With OMP
    !   Time =  2.829 seconds.


! For single presicion
!   
! Without OMP Time =  1.499 seconds.
! 
! With OMP
! Time =  1.406 seconds.

