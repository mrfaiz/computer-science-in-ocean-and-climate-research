program TestMatrix
    implicit none
    real(4),allocatable:: test(:,:),tk(:)
    integer(4):: n,i,j
    n = 3 

    allocate(test(n,n))
    allocate(tk(n))

    call random_number(test)
    
    do i = 1, n
        do j = 1, n 
            print*, test(i,j)
        enddo
    enddo

    tk = test(1,:)

    print*, test(1,:)
    print*,tk

    

    deallocate(test)
    deallocate(tk)
end program TestMatrix
