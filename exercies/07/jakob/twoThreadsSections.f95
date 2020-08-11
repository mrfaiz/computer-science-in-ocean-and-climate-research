program predPreySpace
  use mod_precision
  use omp_lib
  implicit none


  integer(4) :: steps, counter, final_counter
  real(kind=wp) :: t0, T, dt
  real(kind=wp) , allocatable :: preyPop_delta(:), predPop_delta(:), predPop_last(:), preyPop_last(:)
  real(kind=wp), allocatable :: times(:), predPop(:,:), preyPop(:,:), predPop_init(:), preyPop_init(:)
  real(kind=wp), allocatable :: predPop_current(:), preyPop_current(:)
  real(kind=wp), allocatable :: diffPartPred(:), interactPartPred(:), diffPartPrey(:), interactPartPrey(:)
	real(kind=wp), allocatable :: preyPop_intermediate(:), predPop_intermediate(:)
  real(kind=wp) :: alpha, beta, gamma, delta, lambda, mu, epsilon, pi

  logical :: converged, use_improved_euler
  integer(4) :: N, species, save_every_nth, saved_time_steps
  real(kind=wp) :: boxLength, kappa, preFactor, totalLength
  real(kind=wp), allocatable :: diffusionMatrix(:,:)
  integer(4) :: counter1, counter2
  real(8) :: benchStart, benchFinish, differenceNormPred, differenceNormPrey

  namelist /model_parameters/ alpha, beta, gamma, delta, lambda, mu
  namelist /spatial_parameters/ kappa, N
  namelist /time_parameters/ t0, T, dt
  namelist /convergence_parameters/ epsilon
  namelist /export_parameters/ saved_time_steps
	namelist /simulation_parameters/ use_improved_euler

  open (UNIT=10, FILE='predatorprey.nml', STATUS='OLD')
  read (10, NML=model_parameters)
  read (10, NML=spatial_parameters)
  read (10, NML=time_parameters)
  read (10, NML=convergence_parameters)
  read (10, NML=export_parameters)
	read (10, NML=simulation_parameters)

  pi=4.D0*DATAN(1.D0)
  totalLength = 1.0

  steps = int((T-t0)/dt)

  if (steps < int(2*N*N*kappa*(T-t0))) then
    steps = int(2*N*N*kappa*(T-t0))
    write(*,*) "Your step count has been adjusted to ", steps
  end if

  save_every_nth = steps/saved_time_steps

  ! Calculate the temporal divisions
  dt = (T - t0)/DFLOAT(steps)



  ! Calculate the spacial divisions
  boxLength = totalLength/DFLOAT(N)

  ! Prepare initial populations
  allocate(predPop_init(N))
  allocate(preyPop_init(N))

  do counter1 = 1,N
    predPop_init(counter1) = 1.0
    preyPop_init(counter1) = SIN(DFLOAT(counter1)*pi/DFLOAT(N))
  end do

  ! Prepare the diffusionMatrix
  preFactor = kappa/(boxLength*boxLength)
  allocate(diffusionMatrix(N, N))

  counter1 = 1
  do while (counter1 <= N)
    counter2 = 1
    do while (counter2 <= N)
      if (counter1 == counter2) then
        if (counter1 == 1 .OR. counter1 == N) then
          diffusionMatrix(counter1, counter2) = (-1.0)*preFactor
        else
          diffusionMatrix(counter1, counter2) = (-2.0)*preFactor
        end if
      else
        if (counter1 == counter2 + 1 .OR. counter2 == counter1 + 1) then
          diffusionMatrix(counter1, counter2) = preFactor
        else
          diffusionMatrix(counter1, counter2) = 0.0
        end if
      end if
      counter2 = counter2 + 1
    end do
    counter1 = counter1 + 1
  end do


  ! Prepare Arrays for bookkeeping
  allocate(preyPop(steps, N))
  allocate(predPop(steps, N))
  allocate(times(steps))

  preyPop(1, :) = preyPop_init(:)
  predPop(1, :) = predPop_init(:)
  times(1)      = t0

  ! counter = 1
  ! do while (counter <= N)
  !   write(*,*) predPop(1,counter)
  !   counter = counter + 1
  ! end do

  ! Prepare intermediate variables
  allocate(preyPop_delta(N))
  allocate(predPop_delta(N))
  allocate(preyPop_last(N))
  allocate(predPop_last(N))
  allocate(diffPartPred(N))
  allocate(interactPartPred(N))
  allocate(diffPartPrey(N))
  allocate(interactPartPrey(N))
  allocate(predPop_current(N))
  allocate(preyPop_current(N))
	allocate(predPop_intermediate(N))
	allocate(preyPop_intermediate(N))

  predPop_last(:) = predPop_init(:)
  preyPop_last(:) = preyPop_init(:)

	! Initial bookkeeping
	predPop(1,:) = predPop_init(:)
	preyPop(1,:) = preyPop_init(:)
	times(1) = t0 + (DFLOAT(0)*dt)
	final_counter = 1

  ! Main simulation loop
  benchStart = omp_get_wtime()
  counter = 2
  converged = .FALSE.
  do while (counter <= steps .AND. (.NOT. converged))
    !$OMP PARALLEL SECTIONS
      !$OMP SECTION
        preyPop_delta = MATMUL(diffusionMatrix, preyPop_last) + preyPop_last*(alpha - beta*predPop_last - lambda*preyPop_last)
      !$OMP SECTION
        predPop_delta = MATMUL(diffusionMatrix, predPop_last) + predPop_last*(delta*preyPop_last - gamma - mu*predPop_last)
    !$OMP END PARALLEL SECTIONS


		! Improved Euler
		if (use_improved_euler) then
			predPop_intermediate = predPop_last + (dt/2.0)*predPop_delta
			preyPop_intermediate = preyPop_last + (dt/2.0)*preyPop_delta

			!$OMP PARALLEL SECTIONS
	      !$OMP SECTION
	        preyPop_delta = MATMUL(diffusionMatrix, preyPop_intermediate) &
					+ preyPop_intermediate*(alpha - beta*predPop_intermediate - lambda*preyPop_intermediate)
	      !$OMP SECTION
	        predPop_delta = MATMUL(diffusionMatrix, predPop_intermediate) &
					  + predPop_intermediate*(delta*preyPop_intermediate - gamma - mu*predPop_intermediate)
	    !$OMP END PARALLEL SECTIONS
		end if

		! Explicit Euler
		predPop_current = predPop_last + dt*predPop_delta
		preyPop_current = preyPop_last + dt*preyPop_delta

    ! Computing the differences
    differenceNormPred = SQRT(DOT_PRODUCT(predPop_last - predPop_current, predPop_last - predPop_current))
    differenceNormPrey = SQRT(DOT_PRODUCT(preyPop_last - preyPop_current, preyPop_last - preyPop_current))
    converged = (differenceNormPred <= epsilon) .AND. (differenceNormPrey <= epsilon)
    if (converged) then
      write(*,*) "Simulation has converged at T = ", t0 + (DFLOAT(counter)*dt)
      write(*,*) "Step count: ", counter
    end if

    ! Save populations for the next iteration
    predPop_last = predPop_current
    preyPop_last = preyPop_current

    ! Bookkeeping
    if (MOD(counter, save_every_nth) == 0 .AND. (counter/save_every_nth) <= steps) then
       ! write(*,*) "differenceNormPred:", differenceNormPred, "differenceNormPrey:", differenceNormPrey
       predPop(1+(counter/save_every_nth),:) = predPop_current(:)
       preyPop(1+(counter/save_every_nth),:) = preyPop_current(:)
       times(1+(counter/save_every_nth)) = t0 + (DFLOAT(counter)*dt)
       final_counter = counter/save_every_nth
       ! write(*,*) "counter:", counter, "counter/save_every_nth:", (counter/save_every_nth)
     end if

    counter = counter + 1
  end do

  benchFinish = omp_get_wtime()

  open(unit = 20, file = 'results.txt', action = 'write')

  write(20,fmt='(A, ", ")',advance='no') "time"
  counter1 = 1
  do while (counter1 < N)
     write(20,fmt='("pred_",(E18.12), ", " )',advance='no') (DFLOAT(counter1)*boxLength-(boxLength/2.0))
     write(20,fmt='("prey_", (E18.12), ", " )',advance='no') (DFLOAT(counter1)*boxLength-(boxLength/2.0))
     counter1 = counter1 + 1
  end do
  write(20,fmt='("pred_",(E18.12), ", " )',advance='no') (DFLOAT(counter1)*boxLength-(boxLength/2.0))
  write(20,fmt='("prey_", (E18.12), ", " )',advance='no') (DFLOAT(counter1)*boxLength-(boxLength/2.0))
  write(20,*) ""

  counter = 1
  do while (counter <= steps/save_every_nth .AND. counter < final_counter)
    if (MOD(counter, 100) == 0) then
      write(*,*) "Writing: ", ((counter*100)/steps), "%"
    end if
    write(20,fmt='((E25.15), ", ")',advance='no') times(counter)
    counter1 = 1
    do while(counter1 < N)
      write(20,fmt='((E25.15), ", ")',advance='no') predPop(counter,counter1)
      write(20,fmt='((E25.15), ", ")',advance='no') preyPop(counter,counter1)
      counter1 = counter1 + 1
    end do
    write(20,fmt='((E25.15), ", ")',advance='no') predPop(counter,counter1)
    write(20,fmt='((E25.15))',advance='no') preyPop(counter,counter1)
    write(20,*) ""
    counter = counter + 1
  end do

  write(*,*)
  write(*,*) "Benchmark result: ", (benchFinish-benchStart)
  write(*,*)

	deallocate(predPop_intermediate)
	deallocate(preyPop_intermediate)
  deallocate(preyPop_current)
  deallocate(predPop_current)
  deallocate(times)
  deallocate(preyPop_delta)
  deallocate(predPop_delta)
  deallocate(preyPop_last)
  deallocate(predPop_last)
  deallocate(predPop)
  deallocate(preyPop)
  deallocate(diffusionMatrix)
  deallocate(predPop_init)
  deallocate(preyPop_init)
end
