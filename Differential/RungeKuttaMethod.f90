program RungeKuttaMethod
    implicit none
    !dy/dx=y-2x/y; y0=1; x∈(0,1)
    !analytical solution: y=(1+2x)^0.5

    double precision a, b, h, x, k1, k2, k3, k4, y, xk, yk, buff
    double precision :: bb(1:4) = (/1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0/)
    double precision :: cc(1:4) = (/0.0, 0.5, 0.5, 1.0/)
    double precision :: aa(4,4) 
    double precision :: k(4) = (/0,0,0,0/) 
    integer N, i, key, j, m
    
    data aa /0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0/
    a = 0
    b = 1
    y = 1 !initial value
    N = 5
    h = (b-a)/N
    yk = 0
    buff = 0
    read(*,*)key
    write(*,*)' x        y'
    write(*,"(F5.2, F10.6)") a, y
    x = a

    if (key == 1) then
        do i = 1, N
            k1 = y - 2 * x / y
            k2 = y + h*k1/2 - (2*x+h) / (y+h*k1/2)
            k3 = y + h*k2/2 - (2*x+h) / (y+h*k2/2)
            k4 = y + h*k3 - 2*(x+h) / (y+h*k3)
            y = y + h*(k1+2*k2+2*k3+k4)/6
            x = x + h 
            write(*,"(F5.2, F10.6)") x, y
        end do
    else
        do i = 1, N
            do j = 1, 4
                xk = x + h*cc(j)
                do m = 1, j
                    yk = yk + aa(j,m)*k(m)
                end do
                yk = h*yk + y
                k(j) = yk - 2*xk/yk !change the form of f(x,y) to calculate different ODE
                buff = buff + bb(j)*k(j)
                yk = 0
            end do
            x = x + h
            y = y + buff * h
            buff = 0
            write(*,"(F5.2, F10.6)") x, y
        end do
    end if
    stop
end program RungeKuttaMethod
