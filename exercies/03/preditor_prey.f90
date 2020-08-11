program preditor_prey
    !implicit none

    integer(4) :: cycles
    real(8) :: updated_chickens,updated_foxes,model_chicken
    real(8) :: chicken_birth_rate, chicken_death_rate, fox_birth_rate, fox_death_rate, delta_time
    real(8), allocatable :: times(:), foxes(:), chickens(:)

    chicken_birth_rate = 0.5
    chicken_death_rate = 0.015
    fox_birth_rate = 0.015
    fox_death_rate = 0.5

    delta_time = 0.01
    cycles = 4500

    allocate(times(cycles))
    allocate(foxes(cycles))
    allocate(chickens(cycles))

    times(1) = 1
    chickens(1) = 100 
    foxes(1) = 10 

    open(unit = 20, file = 'output.txt', action = 'write')
    !write header
    write(20,fmt='(A, ", ", A, ", ", A)') 'time','foxes','chickens'
    ! Write initala value
    write(20,fmt='(E25.15, ",", E25.15, ",", E25.15)') times(1),foxes(1),chickens(1)
    
    do k = 2, cycles
        times(k) = k
        model_chicken = chicken_birth_rate * chickens(k-1) -  chicken_death_rate * foxes(k-1) * chickens(k-1)
        updated_chickens = chickens(k-1) + delta_time * model_chicken
        updated_foxes = foxes(k-1) + delta_time * ( -fox_death_rate * foxes(k-1) + fox_birth_rate * foxes(k-1) * chickens(k-1))

        write(20,fmt='(E25.15, ",", E25.15, ",", E25.15)') times(k),updated_foxes,updated_chickens

        foxes(k) = updated_foxes
        chickens(k) = updated_chickens
    end do

    close(20)

    deallocate(foxes)
    deallocate(chickens)
    deallocate(times)
   
end program preditor_prey