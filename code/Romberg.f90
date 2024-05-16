program CompositeTrapzoidalRule
    !calculate f(x)=x^2+2x+2 (x:0-1)
    implicit none
    integer :: n, i
    double precision :: a, b, h, Tn, T2n, sum, criteria1, criteria2
    double precision :: arry1(0:6), arry2(0:5)

    a = 0.0
    b = 1.0
    T2n = 0.35  !=(f(a) + f(b)) / 2
    
    sum = 0
    n = 0
    criteria1 = 1
    criteria2 = 1
    arry1(0) = T2n
    
    write(*,"(I1)", advance = 'No') 1
    write(*,"(F15.8)") arry1(0)

    do while ((abs(criteria1) > 0.5D-4 .or. abs(criteria2) > 0.5D-4) .and. n <= 5)
        Tn = T2n
        n = n + 1
        h = (b-a)/2**n
        do i=1, 2**(n-1)
            sum = sum + 1/((a+(2*i-1)*h)**2 + (a+(2*i-1)*h)*2 + 2) !sum = sum + f(a+(2*i-1)*h)
        end do
        T2n = 0.5*Tn + h*sum
        sum = 0

        if (mod(n,2) == 1) then
            arry2(0) = T2n
            do i = 1, n
                arry2(i) = (4**i * arry2(i-1) - arry1(i-1)) / (4**i-1)
            end do
            write(*,"(I1)", advance = 'No') n+1
            do i = 0, n
                write(*,"(F15.8)", advance = 'No') arry2(i)
            end do
            write(*,*)
            criteria1 = arry2(n) - arry1(n-1)
            criteria2 = arry2(n) - arry2(n-1)
        else
            arry1(0) = T2n
            do i = 1, n
                arry1(i) = (4**i * arry1(i-1) - arry2(i-1)) / (4**i-1)
            end do
            write(*,"(I1)", advance = 'No') n+1
            do i = 0, n
                write(*,"(F15.8)", advance = 'No') arry1(i)
            end do
            write(*,*)
            criteria1 = arry1(n) - arry2(n-1)
            criteria2 = arry1(n) - arry1(n-1)
        end if
    end do
    
    stop
end program CompositeTrapzoidalRule