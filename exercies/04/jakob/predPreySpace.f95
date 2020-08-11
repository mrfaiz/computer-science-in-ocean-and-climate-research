program predPreySpace
  !
  implicit none


  integer(4) :: steps, counter
  real(8) :: t_init, t_final, t_delta
  real(8) , allocatable :: preyPop_delta(:), predPop_delta(:), predPop_last(:), preyPop_last(:)
  real(8), allocatable :: times(:), predPop(:,:), preyPop(:,:), predPop_init(:), preyPop_init(:)
  real(8) :: alpha, beta, gamma, delta, lambda, mu

  integer(4) :: spaceDivisions
  real(8) :: boxLength, diffusionCoefficient, preFactor, totalLength
  real(8), allocatable :: diffusionMatrix(:,:)
  integer(4) :: counter1, counter2

  ! Set model and simulation parameters
  spaceDivisions = 100
  totalLength = 1.0
  alpha = 1.0
  beta = 0.1
  lambda = 0.001

  delta = 0.01
  gamma = 0.2
  mu = 0.001
  diffusionCoefficient = 0.001

  t_init = 0.0
  t_final = 100.0
  steps = 10000

  ! Calculate the temporal dimensions
  t_delta = (t_final - t_init)/DFLOAT(steps)

  ! Calculate the spacial dimensions
  boxLength = totalLength/DFLOAT(spaceDivisions)

  ! Prepare initial populations
  allocate(predPop_init(spaceDivisions))
  allocate(preyPop_init(spaceDivisions))


  ! Prepare the initial populations
  counter1 = 1
  do while (counter1 <= spaceDivisions/2)
    predPop_init(counter1) = 0.2
    preyPop_init(counter1) = 1.0
    counter1 = counter1 + 1
  end do

  do while (counter1 <= spaceDivisions)
    predPop_init(counter1) = 2.0
    preyPop_init(counter1) = 0.1
    counter1 = counter1 + 1
  end do

  ! Prepare the diffusionMatrix
  preFactor = diffusionCoefficient/(boxLength*boxLength)
  allocate(diffusionMatrix(spaceDivisions, spaceDivisions))

  counter1 = 1
  do while (counter1 <= spaceDivisions)
    counter2 = 1
    do while (counter2 <= spaceDivisions)
      if (counter1 == counter2) then
        if (counter1 == 1 .OR. counter1 == spaceDivisions) then
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
  allocate(preyPop(steps, spaceDivisions))
  allocate(predPop(steps, spaceDivisions))
  allocate(times(steps))

  preyPop(1, :) = preyPop_init(:)
  predPop(1, :) = predPop_init(:)
  times(1)      = t_init

  ! Prepare intermediate variables
  allocate(preyPop_delta(spaceDivisions))
  allocate(predPop_delta(spaceDivisions))
  allocate(preyPop_last(spaceDivisions))
  allocate(predPop_last(spaceDivisions))

  predPop_last(:) = predPop_init(:)
  preyPop_last(:) = preyPop_init(:)

  ! Main simulation loop
  counter = 2
  do while (counter <= steps)

    ! Calculate the differences
    preyPop_delta = MATMUL(diffusionMatrix, preyPop_last) + preyPop_last*(alpha - beta*predPop_last - lambda*preyPop_last)
    predPop_delta = MATMUL(diffusionMatrix, predPop_last) + predPop_last*(delta*preyPop_last - gamma - mu*predPop_last)

    preyPop_last = preyPop_last + t_delta*preyPop_delta
    predPop_last = predPop_last + t_delta*predPop_delta

    ! Do the bookkeeping
    ! Write current time to array
    times(counter) = t_init + (DFLOAT(counter)*t_delta)
    ! write(*,fmt='(A)',advance='no') "Time: "
    ! write(*,*) times(counter)
    ! write(*,*) counter
    predPop(counter,:) = predPop_last(:)
    preyPop(counter,:) = preyPop_last(:)

    counter = counter + 1
  end do

  ! Write the final results to a file
  open(unit = 20, file = 'results.txt', action = 'write')

  write(20,fmt='(A, ", ")',advance='no') "time"
  counter1 = 1
  do while (counter1 < spaceDivisions)
     write(20,fmt='("pred_",(E18.12), ", " )',advance='no') (DFLOAT(counter1)*boxLength-(boxLength/2.0))
     write(20,fmt='("prey_", (E18.12), ", " )',advance='no') (DFLOAT(counter1)*boxLength-(boxLength/2.0))
     counter1 = counter1 + 1
  end do
  write(20,fmt='("pred_",(E18.12), ", " )',advance='no') (DFLOAT(counter1)*boxLength-(boxLength/2.0))
  write(20,fmt='("prey_", (E18.12), ", " )',advance='no') (DFLOAT(counter1)*boxLength-(boxLength/2.0))
  write(20,*) ""

  counter = 1
  do while (counter <= steps)
    write(20,fmt='((E25.15), ", ")',advance='no') times(counter)
    counter1 = 1
    do while(counter1 < spaceDivisions)
      write(20,fmt='((E25.15), ", ")',advance='no') predPop(counter,counter1)
      write(20,fmt='((E25.15), ", ")',advance='no') preyPop(counter,counter1)
      counter1 = counter1 + 1
    end do
    write(20,fmt='((E25.15), ", ")',advance='no') predPop(counter,counter1)
    write(20,fmt='((E25.15))',advance='no') preyPop(counter,counter1)
    write(20,*) ""
    counter = counter + 1
  end do

  ! deallocate all ressources
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
