real(8), allocatable :: A(:,:), x(:)
real :: res
integer :: n
n = 50
allocate(A(n,n))
call random_number(A)
allocate(x(n))
call random_number(x)

! Trying to access item out of bounds:
print*, x(n+50)
! normal compilation -> prints 0.000000000, compiler doesn't complain
! compiling with -fcheck=bounds leads to a run-time error :
! "Fortran runtime error: Index '3' of dimension 1 of array 'x' above upper bound of 2"

deallocate(x)
deallocate(A)
end
