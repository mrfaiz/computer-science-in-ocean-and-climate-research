program  MatrixMultiplication
    implicit none
    
    real(8), allocatable :: A(:,:), x(:)
    integer :: n 
    n = 10 
    allocate(A(n,n))
    allocate(x(n))

    call random_number(A)
    call random_number(x)

    ! try to access the (n + 1)st element of the vector. What does the compiler say?
    print*, x(11)
    ! compiler showing no error
    ! output is : 0.0000000000000000  

    ! if compile with following command 
    ! "gfortran -fcheck=bound matrixVectorMul.f95"
    ! error found 
    ! "At line 11 of file matrixVectorMul.f95
    ! Fortran runtime error: Index '11' of dimension 1 of array 'x' above upper bound of 10"


    deallocate(A)
    deallocate(x)
    
end program  MatrixMultiplication

