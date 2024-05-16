program CompositeTrapzoidalRule
    implicit none
    integer :: n, i
    double precision :: a, b, h, Tn, T2n, sum

    a = 0.0
    b = 1.0
    Tn = 0
    !T2n = 0.9207355 !=(b-a)(f(a)+f(b))/2
    T2n = 0.35
    sum = 0
    n = 0

    write(*,*) 0.5*T2n+0.5*(sin(0.5))/0.5
    do while (abs(T2n-Tn) > 1.5D-8 .and. n <= 13)
        Tn = T2n
        n = n + 1
        h = (b-a)/2**n
        do i=1, 2**(n-1)
            !sum = sum + (sin(a+(2*i-1)*h))/(a+(2*i-1)*h)
            sum = sum + 1/((a+(2*i-1)*h)**2 + (a+(2*i-1)*h)*2 + 2)
        end do
        T2n = 0.5*Tn + h*sum
        write(*,"('T', I5, '=', F15.10)") 2**n, T2n
        sum = 0
    end do
    write (*,*) n
    
    stop
end program CompositeTrapzoidalRule