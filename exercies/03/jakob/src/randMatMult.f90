program test
  !
  implicit none

  real :: benchCustomn, benchBuiltIn
  real :: x, y, z
  integer :: n

  n = 10

  y = benchCustomn( n )
  z = benchBuiltIn( n )

  write(*,*) 'n -> ', n
  write(*,*) 'benchCustomn( n ) -> ', y
  write(*,*) 'benchBuiltIn( n ) -> ', z

  n = 100

  y = benchCustomn( n )
  z = benchBuiltIn( n )

  write(*,*) 'n -> ', n
  write(*,*) 'benchCustomn( n ) -> ', y
  write(*,*) 'benchBuiltIn( n ) -> ', z

  n = 1000

  y = benchCustomn( n )
  z = benchBuiltIn( n )

  write(*,*) 'n -> ', n
  write(*,*) 'benchCustomn( n ) -> ', y
  write(*,*) 'benchBuiltIn( n ) -> ', z

end

real function benchBuiltIn( n )
  !
  implicit none
  integer :: n
  real(8), allocatable :: A(:,:), x(:), y(:)
  real(8) :: start, finish
  allocate(A(n,n))
  call random_number(A)
  allocate(x(n))
  call random_number(x)

  ! allocate memory for the result vector
  allocate(y(n))
  call cpu_time(start)

  y = MATMUL(A, x)

  call cpu_time(finish)

  deallocate(x)
  deallocate(y)
  deallocate(A)
  benchBuiltIn = finish - start
end function benchBuiltIn

real function benchCustomn( n )
  !
  implicit none
  integer :: n, i, j
  real(8), allocatable :: A(:,:), x(:), y(:)
  real(8) :: start, finish
  allocate(A(n,n))
  call random_number(A)
  allocate(x(n))
  call random_number(x)

  ! allocate memory for the result vector
  allocate(y(n))
  call cpu_time(start)

  do i = 1, n
    y(n) = 0.0
    do j = 1, n
      y(i) = y(i) + x(i)*A(j,i)
    end do
  end do

  call cpu_time(finish)

  deallocate(x)
  deallocate(y)
  deallocate(A)
  benchCustomn = finish - start
end function benchCustomn
