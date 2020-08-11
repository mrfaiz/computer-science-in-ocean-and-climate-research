program mat_multiplication_time_check
    implicit none

    real,allocatable :: A(:,:), x(:)
    integer :: n
    n = 10

    allocate(A(n,n))
    allocate(x(n))
   
    call random_number(A)
    call random_number(x)

    call checkPerformance(A, x, n)
    n = 100 
    call checkPerformance(A, x, n)
    ! n = 1000
    ! call checkPerformance(A, x, n)

    deallocate(A)
    deallocate(x)

contains 

    subroutine checkPerformance(A,x,n)
        implicit none
 
        integer, intent(in) :: n 
        real, intent(in) :: A(n,n), x(n)
        real :: start, finish
        real , allocatable ::  y(:)
        allocate(y(n))

        call cpu_time(start)
        call customeMatrixMul(A,x,n)
        call cpu_time(finish)
        write(*,*) "n =",n,"Time custom multiplication=> ", finish - start
  
        call cpu_time(start)
        y = matmul(A,x)
        call cpu_time(finish)
        write(*,*) "n =",n," Time BuiltIn multiplicatin=> ", finish - start
        deallocate(y)
     end subroutine checkPerformance

    subroutine customeMatrixMul(A, x ,n)
        implicit none
        
        integer, intent(in) :: n
        real, intent(in) :: A(n,n), x(n)
        integer :: i , j
        real , allocatable ::  y(:)
        allocate(y(n))
          do i = 1,n 
            do j = 1,n
                y(i) = y(i) + A(i,j) * x(j) 
            end do
        end do   
        deallocate(y)     
    end subroutine customeMatrixMul

end program mat_multiplication_time_check

! ####################################  Observations  #############################
! $ gfortran  matMulPerformance.f90 
! $ ./a.out 
!  n =          10 Time custom multiplication=>    4.00003046E-06
!  n =          10  Time BuiltIn multiplicatin=>    2.89999880E-05
!  n =         100 Time custom multiplication=>    1.21999998E-04
!  n =         100  Time BuiltIn multiplicatin=>    8.69999640E-05

! for n = 1000 , producing error
! Program received signal SIGSEGV: Segmentation fault - invalid memory reference.

! Backtrace for this error:
! #0  0x7faecaff8cea
! #1  0x7faecaff7e75
! #2  0x7faecae2e46f
! #3  0x559053c934c0
! #4  0x559053c9328e
! #5  0x7faecae0f1e2
! #6  0x559053c9333d
! #7  0xffffffffffffffff
! Segmentation fault (core dumped)

! $ gfortran -O  matMulPerformance.f90 
! $ ./a.out 
!  n =          10 Time custom multiplication=>    2.00001523E-06
!  n =          10  Time BuiltIn multiplicatin=>    1.00000761E-06
!  n =         100 Time custom multiplication=>    3.59999249E-05
!  n =         100  Time BuiltIn multiplicatin=>    1.69998966E-05

! fortran -O1  matMulPerformance.f90 
! $ ./a.out 
!  n =          10 Time custom multiplication=>    2.00001523E-06
!  n =          10  Time BuiltIn multiplicatin=>    1.00000761E-06
!  n =         100 Time custom multiplication=>    3.30000184E-05
!  n =         100  Time BuiltIn multiplicatin=>    2.00000359E-05

! $ gfortran -O2  matMulPerformance.f90 
! $ ./a.out
!  n =          10 Time custom multiplication=>    3.00002284E-06
!  n =          10  Time BuiltIn multiplicatin=>    1.00000761E-06
!  n =         100 Time custom multiplication=>    3.40000261E-05
!  n =         100  Time BuiltIn multiplicatin=>    1.90000283E-05