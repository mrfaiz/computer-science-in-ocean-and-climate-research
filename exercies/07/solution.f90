program PreditorPreySpace
    use precision
    USE OMP_LIB
    implicit none

    integer(4) :: N, calculated_steps,i,j,total_steps
    real :: start_time, end_time
    real(kind=wp) :: alpha,beta,gamma,delta,lambda,mu,kappa,t0,T,dt,delta_single_box_length,diffusion_factor,PI,epsilon,preditor_diff,prey_diff
    
    real(kind=wp), allocatable :: boxes(:),times(:),D(:,:),preditor_matrix(:,:), prey_matrix(:,:)
    real(kind=wp), allocatable :: xk(:), yk(:),prey_f1(:),preditor_f2(:)

    namelist /model_parameters/ alpha, beta, gamma, delta, lambda, mu
    namelist /spatial_parameters/ kappa, N
    namelist /time_parameters/ t0, T, dt

    open(unit = 20, file = 'predatorprey.nml', action = 'read')
    read(20,nml=model_parameters)
    read(20,nml=spatial_parameters)
    read(20,nml=time_parameters)
    close(20)

    ! Constants
    PI=4.D0*DATAN(1.D0) ! from stack overflow
    print*,PI 
    epsilon = 0.00

  
    ! Make sure in your code that the condition n ≥ 2N 2 κT is always satisfied
    calculated_steps = int(2*N*N*kappa*T)

    ! Allocation
    allocate(boxes(N))
    allocate(times(calculated_steps))
    allocate(D(N,N))
    allocate(preditor_matrix(n,N))
    allocate(prey_matrix(calculated_steps,N))
    allocate(yk(N))
    allocate(xk(N))
    allocate(preditor_f2(N))
    allocate(prey_f1(N))
    ! Box Data
    delta_single_box_length = 1.0/DFLOAT(N) ! also h
    do i = 1, N
        boxes(i) = i * delta_single_box_length
    enddo
    
    ! Initial values
    times(1) = 0.0 
    ! populate initial(random data) data inside preditor, prey and times matrix

    do i = 1, N 
        prey_matrix(1,i) = SIN((DFLOAT(i)*PI)/DFLOAT(N))
        preditor_matrix(1,i) = 1.0    
    enddo

    ! Calculate Diffusion multiplication Factor
    diffusion_factor = kappa/(delta_single_box_length * delta_single_box_length)
   
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

    !$OMP DO  
    do i = 2, calculated_steps 
        xk = prey_matrix(i-1,:)
        yk = preditor_matrix(i-1,:)


        
        prey_f1 = matmul(D,xk) + xk * (alpha - (beta * yk) - (lambda * xk))
        preditor_f2 = matmul(D,yk) + yk * ((delta * xk) - gamma - (mu * yk))  

        prey_matrix(i,:) = xk + dt * prey_f1
        preditor_matrix(i,:) = yk + dt * preditor_f2

        times(i) = times(i-1) + dt
    enddo
    !$OMP END DO

    ! contains 

    !     subroutine euclidean_norm(i)
    !         implicit none
    !         integer, intent(in) :: i

    !     end subroutine euclidean_norm

   
end program PreditorPreySpace

