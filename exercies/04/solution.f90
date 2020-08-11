program PreditorPreySpace
    implicit none
    real(8) :: alpha, beta, lamda, delta, gamma, mu,kappa,delta_single_box_length,delta_time,diffusion_factor
    integer(4) :: T, cap_n,i,n,j
    real(8), allocatable :: boxes(:),times(:),D(:,:),preditor_matrix(:,:), prey_matrix(:,:)
    real(8), allocatable :: xk(:), yk(:),prey_f1(:),preditor_f2(:)

    alpha = 1.0
    beta = 0.1

    delta = 0.1
    gamma = 0.2

    lamda = 0.001
    mu = 0.001

    kappa = 0.001
    T = 500
    cap_n = 100
    
    ! Calculate n
    n = 4 * cap_n * cap_n * kappa * T
    print*,n
    ! Allocation
    allocate(boxes(cap_n))
    allocate(times(n))
    allocate(D(cap_n,cap_n))
    allocate(preditor_matrix(n,cap_n))
    allocate(prey_matrix(n,cap_n))
    allocate(yk(cap_n))
    allocate(xk(cap_n))
    allocate(preditor_f2(cap_n))
    allocate(prey_f1(cap_n))
    ! Box Data
    delta_single_box_length = 1.0/DFLOAT(cap_n) ! also h
    do i = 1, cap_n
        boxes(i) = i * delta_single_box_length
    enddo
    
    ! Calculate Delta_t and populate data to time vector
    delta_time = (T - 0.0)/ DFLOAT(n)
    
    ! Calculate Diffusion multiplication Factor
    diffusion_factor = kappa/(delta_single_box_length * delta_single_box_length)
    ! print*, diffusion_factor
   
    ! populate data in diffussion matrix
    do i = 1, cap_n ! row
        do j = 1 , cap_n !colum
            if (i==j) then
                if (i == 1 .OR. i == cap_n) then
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

  

    do i = 1, cap_n 
        if(i<=(cap_n/2)) then
            preditor_matrix(1,i) = 0.2
            prey_matrix(1,i) = 1.0
        else
            preditor_matrix(1,i) = 2.0
            prey_matrix(1,i) = 0.1
        endif
    enddo

    times(1) = 0.0 

    ! Main Loop
    do i = 2, n 
        xk = prey_matrix(i-1,:)
        yk = preditor_matrix(i-1,:)

        ! f 1 (x k , y k ) = Dx k + x k ∗ (α − βy k − λx k )
        ! f 2 (x k , y k ) = Dy k + y k ∗ (δx k − γ − μy k )
        prey_f1 = matmul(D,xk) + xk * (alpha - (beta * yk) - (lamda * xk))
        preditor_f2 = matmul(D,yk) + yk * ((delta * xk) - gamma - (mu * yk))  

        prey_matrix(i,:) = xk + delta_time * prey_f1
        preditor_matrix(i,:) = yk + delta_time * preditor_f2

        ! tk+1 = tk + delta_t * k
        times(i) = times(i-1) + delta_time
    enddo

    ! Write rest of the list
    open(unit = 20, file = '/home/faiz/SS_2020/Ocean/exercies/04/prey.txt', action = 'write')
    open(unit = 21, file = '/home/faiz/SS_2020/Ocean/exercies/04/preditor.txt', action = 'write')
    open(unit = 22, file = '/home/faiz/SS_2020/Ocean/exercies/04/time.txt', action = 'write')

    do i = 1,n 
        do j = 1, cap_n 
            write(20,fmt='(E25.15,",")',advance='no') prey_matrix(i,j)
            write(21,fmt='(E25.15,",")',advance='no') preditor_matrix(i,j)
        enddo
        write(20, *) ''  ! this gives you the line break
        write(21, *) ''  ! this gives you the line break
        write(22,fmt='(E25.15)') times(i)
    enddo

 
    close(20)
    close(21)
    close(22)

    !Deallocation
    deallocate(boxes)
    deallocate(times) 
    deallocate(D)
    deallocate(preditor_matrix)
    deallocate(prey_matrix)
    deallocate(yk)
    deallocate(xk)
    deallocate(preditor_f2)
    deallocate(prey_f1)
end program PreditorPreySpace