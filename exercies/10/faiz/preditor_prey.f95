program PreditorPreySpace
    use precision
    USE OMP_LIB
    implicit none

    integer(4) :: N, steps, calculate_steps,counter,previous_index

    real(kind=wp) :: alpha,beta,gamma,delta,lambda,mu,kappa,t0,T,dt
    real(kind=wp) :: xk_preditor,yk_prey,prey_delta_t,prey_intermediate_delta_t,preditor_delta_t, preditor_intermediate_delta_t
    real(kind=wp) ::  xk_plus_one_preditor, yk_plus_one_prey
    
    real(kind=wp), allocatable :: times(:),preditor(:), prey(:)
    logical :: use_improved_euler

    namelist /model_parameters/ alpha, beta, gamma, delta, lambda, mu
    namelist /time_parameters/ t0, T, dt
    namelist /simulation_parameters/ use_improved_euler

    character(len = 35)::caption_fmt,file_write_fmt

    ! call cp(start_time)
    open(unit = 20, file = 'predatorprey.nml', action = 'read')
    read(20,nml=model_parameters)
    read(20,nml=time_parameters)
    read(20,nml=simulation_parameters)

    steps = int((T-t0)/dt)

    ! print*,alpha,beta,gamma,delta,lambda,mu 
    
  ! Allocation
    allocate(times(steps))
    allocate(preditor(steps))
    allocate(prey(steps))

    times(1) = t0 
    preditor(1) = 1.0
    prey(1) = 0.2

    caption_fmt = '(A, ", ", A, ", ", A)'
    file_write_fmt = '(E25.15, ",", E25.15, ",", E25.15)'

    open(unit = 20, file = 'output.txt', action = 'write')
    !write header
    write(20,fmt= caption_fmt) 'time','prey','preditor'
    ! Write initala value
    write(20,fmt=file_write_fmt) times(1),prey(1),preditor(1)
    print*, use_improved_euler   
    do counter = 2 , steps
        previous_index = counter-1

        times(counter) = times(previous_index) + dt
        
        xk_preditor = preditor(previous_index)
        yk_prey = prey(previous_index)

        

        preditor_delta_t = xk_preditor * (alpha - beta*yk_prey - lambda*xk_preditor)
        prey_delta_t = yk_prey * (delta*xk_preditor - gamma - mu*yk_prey)

        ! imporved euler start 
        if (use_improved_euler) then

          preditor_intermediate_delta_t = xk_preditor + (dt/2) * preditor_delta_t
          prey_intermediate_delta_t = yk_prey + (dt/2) * prey_delta_t

          preditor_delta_t =  preditor_intermediate_delta_t * (alpha - (beta*prey_intermediate_delta_t) - lambda*&
          preditor_intermediate_delta_t)
          prey_delta_t = prey_intermediate_delta_t * (delta*preditor_intermediate_delta_t - gamma - mu*prey_intermediate_delta_t)

         ! imporved euler end
        endif

        xk_plus_one_preditor = xk_preditor + dt * preditor_delta_t
        yk_plus_one_prey  = yk_prey + dt * prey_delta_t

        write(20,fmt=file_write_fmt) times(counter),yk_plus_one_prey,xk_plus_one_preditor
        preditor(counter) = xk_plus_one_preditor
        prey(counter) = yk_plus_one_prey
        ! write(*,*) xk_plus_one_preditor,yk_plus_one_prey

        ! for improved euler
        ! preditor_intermediate_delta_t = 


    enddo




    ! DeAllocation
    deallocate(times)
    deallocate(preditor)
    deallocate(prey)
end program PreditorPreySpace

