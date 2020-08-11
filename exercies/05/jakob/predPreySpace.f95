program predPreySpace
  use mod_precision
  use omp_lib
  implicit none


  integer(4) :: steps, counter
  real(kind=wp) :: t0, T, dt
  real(kind=wp) , allocatable :: preyPop_delta(:), predPop_delta(:), predPop_last(:), preyPop_last(:)
  real(kind=wp), allocatable :: times(:), predPop(:,:), preyPop(:,:), predPop_init(:), preyPop_init(:)
  real(kind=wp) :: alpha, beta, gamma, delta, lambda, mu

  integer(4) :: N, species
  real(kind=wp) :: boxLength, kappa, preFactor, totalLength
  real(kind=wp), allocatable :: diffusionMatrix(:,:)
  integer(4) :: counter1, counter2
  real(8) :: benchStart, benchFinish

  namelist /model_parameters/ alpha, beta, gamma, delta, lambda, mu
  namelist /spatial_parameters/ kappa, N
  namelist /time_parameters/ t0, T, dt

  open (UNIT=10, FILE='predatorprey.nml', STATUS='OLD')
  read (10, NML=model_parameters)
  read (10, NML=spatial_parameters)
  read (10, NML=time_parameters)

  totalLength = 1.0

  steps = int((T-t0)/dt)

  if (steps < int(2*N*N*kappa*T)) then
    steps = int(2*N*N*kappa*T)
    write(*,*) "Your step count has been adjusted"
  end if

  ! Calculate the temporal divisions
  dt = (T - t0)/DFLOAT(steps)

  ! Calculate the spacial divisions
  boxLength = totalLength/DFLOAT(N)

  ! Prepare initial populations
  allocate(predPop_init(N))
  allocate(preyPop_init(N))


  counter1 = 1
  do while (counter1 <= N/2)
    predPop_init(counter1) = 0.2
    preyPop_init(counter1) = 1.0
    counter1 = counter1 + 1
  end do

  do while (counter1 <= N)
    predPop_init(counter1) = 2.0
    preyPop_init(counter1) = 0.1
    counter1 = counter1 + 1
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

  predPop_last(:) = predPop_init(:)
  preyPop_last(:) = preyPop_init(:)

  ! Main simulation loop
  benchStart = omp_get_wtime()
  counter = 2
  species = 0
!$OMP PARALLEL PRIVATE( species ) num_threads(2)
  do counter = 2, steps
    ! Calculate the differences
    !if (MOD(counter, 100) == 0) then
      !write(*,*) "Simulation: ", ((counter*100)/steps), "%"
    !end if

    !$OMP DO
    do species = 0, 1
      if (species == 0) then
        preyPop_delta = MATMUL(diffusionMatrix, preyPop_last) + preyPop_last*(alpha - beta*predPop_last - lambda*preyPop_last)
        preyPop_last = preyPop_last + dt*preyPop_delta
        preyPop(counter,:) = preyPop_last(:)
      else
        predPop_delta = MATMUL(diffusionMatrix, predPop_last) + predPop_last*(delta*preyPop_last - gamma - mu*predPop_last)
        predPop_last = predPop_last + dt*predPop_delta
        predPop(counter,:) = predPop_last(:)
      end if
    end do
    !$OMP END DO

    times(counter) = t0 + (DFLOAT(counter)*dt)

    !counter = counter + 1
  end do
!$OMP END PARALLEL

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
  do while (counter <= steps)
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
