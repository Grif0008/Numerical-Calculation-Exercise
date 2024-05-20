program RungeKuttaMethod
    implicit none
    !y'=-y+cos2x-2sin2x+2xe^-x; y0=1; x∈(0,2)
    !analytical solution: (x^2)*exp(-x)+cos(2*x)

    double precision a, b, h, x, y, k1, k2, k3, k4
    double precision x_rk(4), y_rk(4), fn(4)
    double precision p(2), c(2)
    double precision m, fm
    real, parameter :: e = 2.718281828459045
    integer N, i, j
    integer :: fileid = 100
    character(len = 10) :: filename = 'data.txt' 
    
    a = 0
    b = 2
    y = 1 !initial value
    p(1) = 0.0
    c(1) = 0.0
    write(*,*) 'Please Enter an Integer Number for N (bigger than 4):'
    read (*,*) N
    h = (b-a)/N
  
    open(unit = fileid, file = filename)

    write(fileid,*)'   x        y'
    write(fileid,"(F8.4, F10.6)") a, y
    x = a

    !Use normal Runge-Kutta Method to calculate the first 3 values.
    y_rk(1) = y
    x_rk(1) = x
    do i = 1, 3
        k1 = -y + cos(2*x) - 2*sin(2*x) + 2*x*(e**(-x))
        k2 = y + h*k1/2 - (2*x+h) / (y+h*k1/2)
        k3 = y + h*k2/2 - (2*x+h) / (y+h*k2/2)
        k4 = y + h*k3 - 2*(x+h) / (y+h*k3)
        y = y + h*(k1+2*k2+2*k3+k4)/6
        y_rk(i+1) = y
        x = x + h 
        x_rk(i+1) = x
        write(fileid,"(F8.4, F10.6)") x, y
    end do

    !Use Adams-Predictior-Corrector Method to finish the following calculation.
    do i = 4, N
        do j = 1, 4
            fn(j) = -y_rk(j) + cos(2*x_rk(j)) - 2*sin(2*x_rk(j)) + 2*x_rk(j)*(e**(-x_rk(j)))
        end do
        
        p(2) = y + h * (55*fn(4) - 59*fn(3) + 37*fn(2) - 9*fn(1)) / 24

        m = p(2) + 251 * (c(1) - p(1)) / 270
        x = x + h
        fm = -m + cos(2*x) - 2*sin(2*x) + 2*x*(e**(-x))

        c(2) = y + h * (9*fm + 19*fn(4) - 5*fn(3) + fn(2)) / 24
        y = c(2) - 19 * (c(2) - p(2)) / 270
        
        do j = 1, 3
            y_rk(j) = y_rk(j+1)
            x_rk(j) = x_rk(j+1)
        end do
        y_rk(4) = y 
        x_rk(4) = x 
        p(1) = P(2)
        c(1) = c(2)
        write(fileid,"(F8.4, F10.6)") x, y
    end do
    close(fileid)

    write(*,*) 'Done! Check the Results in "data.txt".'

    stop
end program RungeKuttaMethod
