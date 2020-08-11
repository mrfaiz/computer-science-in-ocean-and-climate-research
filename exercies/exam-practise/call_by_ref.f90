program call_by_ref
    implicit none
    integer:: x,y 
    x = 4
    y = 5 

    write(*,*) "before :x = ",x," y= ",y
    call swap(x,y)
    write(*,*) "after :x = ",x," y= ",y
end program call_by_ref

subroutine swap(x,y)
    integer:: x,y,help
    help = x 
    x = y 
    y = help
end subroutine swap
