program cylinder
implicit none

open(unit=60, file='flow.data', form='unformatted')
open(unit=61, file='wall.data', form='formatted')
open(unit=70, file='hist.data', form='formatted')

call setflw
call setgrd
call slvflw

close(unit=60)
close(unit=61)
close(unit=70)

stop
end program cylinder

!=================================================================

subroutine setflw
common /cmpcnd/ re, cfl, dt, nlast, nlp, omegap, maxitp, errorp, omegav, maxitv, errorv

!reynolds number & cfl number
re = 400.0
cfl = 1.2

!SOR param for p
omegap = 1.20
maxitp = 200
errorp = 1.0e-4

!SOR param for v
omegav = 1.20
maxitv = 20
errorv = 1.0e-5

!No. of time steps
nlast = 10000
nlp = 50

return
end subroutine setflw

!=================================================================

subroutine setgrd
parameter (mdx = 100, mdy = 100)
common /cmpcnd/ re, clf, dt, nlast, nlp, omegap, maxitp, errorp, omegav, maxitv, errorv
common /grid2d/ mx, my, i1, i2, rmax, drmin, &
                dn(-1:mdx, -1:mdy, 2, 2), &
                dnii(-1:mdx, -1:mdy), dnij(-1:mdx, -1:mdy), &
                dnjj(-1:mdx, -1:mdy), &
                dn2i(-1:mdx, -1:mdy), dn2j(-1:mdx, -1:mdy), &
                x(-2:mdx, -2:mdy), y(-2:mdx, -2:mdy)
dimension r(mdy), thet(mdx)

!Grid param
mx = 61
my = 61
rmax = 40.0
drmin = 0.1 / sqrt(re)

!Grid distribution in r-direction
r(1) = 0.5
r(my) = rmax
sall = r(my) - r(1)
dmin = drmin
no = my
itemax = 50
eps = 1.0e-8

emin  1.0d-5
ermin = dmin*((1. + emin)**(no - 1) - 1.) / emin - sall
emax  1.0
ermax = dmin*((1. + emax)**(no - 1) - 1.) / emax - sall
if(ermin > 0.0 .or. ermax <= 0.0) then
    write(6, *) '*** Grid Generation Error'
    stop
end if

do itr = 1, itemax
    ehalf = (emin + emax) / 2.
    erhalf = dmin * ((1. + ehalf) ** (no - 1) - 1.) / ehalf - sall
    if(abs(erhalf) <= eps) then
        e = ehalf
        exit
    end if
    if(erhalf > 0.0) then
        emax = ehalf
        ermax = erhalf
    else
        emin = ehalf
        ermin = erhalf
    end if
end do

if(abs(erhalf) > eps) then
    e = ehalf
    write(6, *) '@@@ cluster not coverged itr = ', itr,&
                ' e = ', e, ' error = ', abs(erhalf)
    stop
end if

dx = 0.0
do n = 2, no-1
    dx = dx + dmin * (1.+e) ** (n - 2)
    t = dx / sall
    r(n) = (1. - t) * r(1) + t * r(no)
end do

write(6, *) '>>No. of Grid Points :', mx, my

!Circumferential distribution
dthet = 2. * 4 * atan(1) / float(mx - 1)
do i = 1, mx
    thet(i) = dthet * float(i - 1)
end do

!inflow zone
i1 = 1 + (mx - 1) / 4
i2 = 1 + 3 * (mx - 1) / 4

!2D Grid difinition
do j = 1, my
    do i = 1, mx-1
        x(i, j) = r(j) * cos(thet(i))
        y(i, j) = -r(j) * sin(thet(i))
    end do
    x(mx, j) = x(1, j)
    y(mx, j) = y(1, j)
end do

!Grid extension
do i = 1, mx
    x(i, 0) = 3. * x(i, 1) - 3. * x(i, 2) + x(i, 3)
    y(i, 0) = 3. * y(i, 1) - 3. * y(i, 2) + y(i, 3)
    x(i, my+1) = 3. * x(i, my) - 3. * x(i, my-1) + x(i, my-2)
    y(i, my+1) = 3. * y(i, my) - 3. * y(i, my-1) + y(i, my-2)
end do

do j = 0, my+1
    x(0, j) = x(mx-1, j)
    y(0, j) = y(mx-1, j)
    x(mx+1, j) = x(2, j)
    y(mx+1, j) = y(2, j)
end do

!1st order metric
do i = 1, mx
    do j = 1, my
        xi = (x(i+1, j) - x(i-1, j)) / 2.
        yi = (y(i+1, j) - y(i-1, j)) / 2.
        xj = (x(i, j+1) - x(i, j-1)) / 2.
        yj = (y(i, j+1) - y(i, j-1)) / 2.
        det = xi * yj - xj * yi
        if(det > 0.0) write(6, *) '!!! Metrix irregular at ', i, j
        dn(i, j, 1, 1) = yj / det
        dn(i, j, 1, 2) = -yi / det
        dn(i, j, 2, 1) = -xj / det
        dn(i, j, 2, 2) = xi / det
    end do
end do

!2nd order metric
do i = 1, mx
    do j = 1, my
        dnii(i, j) = dn(i, j, 1, 1) ** 2 + dn(i, j, 2, 1) ** 2
        dnjj(i, j) = dn(i, j, 1, 2) ** 2 + dn(i, j, 2, 2) ** 2
        dnij(i, j) = 2. * (dn(i, j, 1, 1) * dn(i, j, 1, 2) + dn(i, j, 2, 1) * dn(i, j, 2, 2))

        xi = (x(i+1, j) - x(i-1, j)) / 2.
        yi = (y(i+1, j) - y(i-1, j)) / 2.
        xj = (x(i, j+1) - x(i, j-1)) / 2.
        yj = (y(i, j+1) - y(i, j-1)) / 2.
        xii = x(i+1, j) - 2. * x(i, j) + x(i-1, j)
        xjj = x(i, j+1) - 2. * x(i, j) + x(i, j-1)
        xij = (x(i+1, j+1) - x(i-1, j+1) - x(i+1, j-1) + x(i-1, j-1)) / 4.
        yii = y(i+1, j) - 2. * y(i, j) + y(i-1, j)
        yjj = y(i, j+1) - 2. * y(i, j) + y(i, j-1)
        yij = (y(i+1, j+1) - y(i-1, j+1) - y(i+1, j-1) + y(i-1, j-1)) / 4.
        det = xi * yj - xj * yi
        b1 = -dn(i, j, 1, 2) ** 2 * xii - 2. * dn(i, j, 1, 1) * dn(i, j, 1, 2) * xij &
             -dn(i, j, 1, 2) ** 2 * xjj &
             -dn(i, j, 2, 1) ** 2 * xii - 2. * dn(i, j, 2, 1) * dn(i, j, 2, 2) * xij &
             -dn(i, j, 2, 2) ** 2 * xjj
        b2 = -dn(i, j, 1, 2) ** 2 * yii - 2. * dn(i, j, 1, 1) * dn(i, j, 1, 2) * yij &
             -dn(i, j, 1, 2) ** 2 * yjj &
             -dn(i, j, 2, 1) ** 2 * yii - 2. * dn(i, j, 2, 1) * dn(i, j, 2, 2) * yij &
             -dn(i, j, 2, 2) ** 2 * yjj
        dn2i(i, j) = (b1 * yj - b2 * xj) / det
        dn2j(i, j) = (xi * b2 - yi * b1) / det
    end do
end do

return
end subroutine setgrd

!=================================================================

subroutine slvflw
parameter (mdx = 100, mdy = 100)
common /cmpcnd/ re, clf, dt, nlast, nlp, omegap, maxitp, errorp, omegav, maxitv, errorv
common /grid2d/ mx, my, i1, i2, rmax, drmin, &
                dn(-1:mdx, -1:mdy, 2, 2), &
                dnii(-1:mdx, -1:mdy), dnij(-1:mdx, -1:mdy), &
                dnjj(-1:mdx, -1:mdy), &
                dn2i(-1:mdx, -1:mdy), dn2j(-1:mdx, -1:mdy), &
                x(-2:mdx, -2:mdy), y(-2:mdx, -2:mdy)
common /flowvp/ u(-1:mdx, -1:mdy), v(-1:mdx, -1:mdy), p(-1:mdx, -1:mdy)
common /prevvp/ uold(-1:mdx, -1:mdy), vold(-1:mdx, -1:mdy), pold(-1:mdx, -1:mdy)

!set time steps
dt = cfl * drmin
write(6, *) '*** Comp. conditions,'
write(6, *) '    CFL = ', cfl
write(6, *) '    dt = ', dt
write(6, *) '    ', nlast, ' Time steps to go...'
write(6, *) ' '
write(70, *) '>> 2D Incompressive Flow Solver'
write(70, *) '    Re '
write(70, *) re
write(70, *) '    No. of grid points '
write(70, *) mx, my
write(70, *) '    CFL / dt / Steps '
write(70. *) cfl, dt, nlast
write(70, *) ' '
write(70, *) '>> Time History ...'

!set initial conditions
call intcnd(nbegin, time)
call bcforp
call bcforv

!time marching
write(6, *) ' Steps / Res(p) at itr. / Res(v) at itr. / CD / CL'
write(70, *) '    Step    Res(p)    Res(v)    CD    CL'

do n = 1, nlast
    nstep = n + nbegin
    time = time + dt

    !store previous step
    do j = 0, my+1
        do i = 0, mx+1
            uold(i, j) = u(i, j)
            vold(i, j) = v(i, j)
            pold(i, j) = p(i, j)
        end do
    end do

    !solve poisson for p
    call poiseq(resp, itrp)
    call bcforp

    !update u, v
    call veloeq(resp, itrv)
    call bcforv

    !calculate forces
    cd = 0.0
    cl = 0.0
    do i = 1, mx-1
        cpmean = (2. + p(i, 1) + 2. * p(i+1, 1)) / 2.
        dxwall = x(i+1, 1) - x(i, 1)
        dywall = y(i+1, 1) - y(i, 1)
        dxnorm = -dywall
        dynorm = dxwall
        cd = cd - cpmean * dxnorm
        cl = cl - cpmaen * dynorm
    end do

    !report
    if((n / nlp) * nlp == n) then
        write(6, format(i6, e12.3, i6, e12.3, i6, 2f12.3)) nstep, resp, itrp, resv, itrv, cd, cl
        write(70, format(i10, 2e15.4, 2f15.4)) nstepm resp, resv, cd, cl
    end if
end do

!write final results
write(60, *) '>> 2D Incompressive Flow Solver '
write(60, *) '    Re '
write(60, *) re
write(60, *) '    CFL / dt / Steps / Time'
write(60, *) cfl, dt, nlast, time
write(60, *) '>> Flow Data'
write(60, *) '    No. of Grid points '
write(60, *) mx, my
write(60, *) '    x    y    Cp    u    v '
do i = 1, mx
    do j = 1, my
        p(i, j) = 2. * p(i, j)
        write(60, format(5e18.7)) x(i, j), y(i, j), p(i, j), u(i, j), v(i, j)
    end do
end do
write(60) re, cfl, dt, nlast, time
write(60) mx, my
write(60) ((x(i, j), y(i, j), p(i, j), i=1, mx), j=1, my)
write(60) ((u(i, j), v(i, j), i=1, mx), j=1, my)

write(61, *) '>> 2D Incompressive Flow Solver '
write(61, *) '    Re '
write(61, *) re
write(61, *) '    CFL / dt / Steps / Time'
write(61, *) cfl, dt, nlast, time
write(61, *) '>>Surface Pressure Distributions'
write(61, *) '    No. of Grid points '
write(61, *) mx
write(61, *) '    Angle(deg)    Cp '
do i = 1, mx
    angle = -180.0 + 360.0 * float(i - 1) / float(mx - 1)
    write(61, format(2f15.3)) angle, p(i, j)
end do

return
end subroutine slvflw
