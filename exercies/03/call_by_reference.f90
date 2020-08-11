program CallByRef
    implicit none
    integer :: x, y
    x = 1
    y = 2
    write(*,*) "before swap: x =>", x , " y = ",y
    call swap(x,y)
    write(*,*) "before swap: x =>", x , " y = ",y
end program CallByRef

subroutine swap(x,y)
    integer:: x,y,help
    help = x
    x = y 
    y = help
end subroutine swap