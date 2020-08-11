program test
  !
  implicit none
  integer(4) :: steps, counter
  real(8) :: t_init, t_final, t_delta, predPop_init, preyPop_init
  real(8) :: preyPop_delta, predPop_delta, predPop_last, preyPop_last
  real(8), allocatable :: times(:), predPop(:), preyPop(:)
  real(8) :: alpha, beta, gamma, delta, lambda, mu

  alpha = 1.0
  beta = 0.1
  lambda = 0.001

  delta = 0.01
  gamma = 0.2
  mu = 0.001

  t_init = 0.0
  t_final = 100.0
  steps = 100000

  predPop_init = 1.0
  preyPop_init = 1.0

  ! Calculate the time delta
  t_delta = (t_final - t_init)/(DFLOAT(steps))

  ! Allocate memory for the results
  allocate(times(steps))
  allocate(predPop(steps))
  allocate(preyPop(steps))

  times(1) = t_init
  predPop(1) = predPop_init
  preyPop(1) = preyPop_init

  predPop_last = predPop_init
  preyPop_last = preyPop_init
  counter = 2

  do while(counter <= steps)
      times(counter) = t_init + (DFLOAT(counter)*t_delta)
      preyPop_delta = preyPop_last*(alpha - beta*predPop_last - lambda*preyPop_last)
      predPop_delta = predPop_last*(delta*preyPop_last - gamma - mu*predPop_last)

      preyPop_last = preyPop_last + t_delta*preyPop_delta
      predPop_last = predPop_last + t_delta*predPop_delta

      predPop(counter) = predPop_last
      preyPop(counter) = preyPop_last

      counter = counter + 1
  end do


  open(unit = 20, file = 'predPreyResults.txt', action = 'write')
  !write header
  write(20,fmt='(A, ", ", A, ", ", A)') 'time','predatorPopulation','preyPopulation'


  !write data
  do counter=1, steps
      write(20,fmt='(E25.15, ",", E25.15, ",", E25.15)') times(counter),predPop(counter),preyPop(counter)
  end do

  close(20)

  deallocate(times)
  deallocate(predPop)
  deallocate(preyPop)
end
