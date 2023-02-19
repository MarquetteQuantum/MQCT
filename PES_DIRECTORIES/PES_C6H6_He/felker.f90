module params
        implicit none
        real*8, parameter :: r0 = 7.2205270d0, a = 0.5724460d0
        real*8, parameter :: W0 = 0.0237670d0, bz = 1.2648900d0
        real*8, parameter :: c3 = 0.0583870d0, c4 = -6.9148510d0
        real*8, parameter :: c5 = -2.0838080d0, c6 = 85.4295680d0
        real*8, parameter :: c11 = -16.9569970d0, c12 = 0.5864290d0
        real*8, parameter :: c22 = 9.0731840d0, c13 = 1.8219200d0
        real*8, parameter :: c14 = 1.8747760d0, c23 = -0.0983960d0
        real*8, parameter :: c111 = -3.2351220d0, c112 = -2.6315040d0
        real*8, parameter :: c122 = 1.4237770d0, c113 = -2.1171520d0
        real*8 c0 
        real*8 x(6), y(6), z(6)	
		real*8, parameter :: shift_bikram = 212.969407615856910d0
end module

subroutine init
        use params
        implicit none
        c0 = 6.0d0*(1.0d0 + c3 + c4 + c5 + c6) &
                + 15.0d0*(c11 + c22 + 2.0d0*(c12 + c13 + c14 + c23)) &
                + 20.0d0*(c111 + 3.0d0*(c112 + c113 + c122))
        x(1) =   1.203775311260370d0
        x(2) =   1.203775311260370d0
        x(3) =   0.00d0
        x(4) =  -1.203775311260370d0
        x(5) =  -1.203775311260370d0
        x(6) =   0.00d0

        y(1) =  -0.695d0
        y(2) =   0.695d0
        y(3) =   1.390d0
        y(4) =   0.695d0
        y(5) =  -0.695d0
        y(6) =  -1.390d0

        z(:) =   0.00d0
end subroutine

subroutine felker_pes(R, beta, gamma, V)
        use params
        implicit none
        integer k, l, m
        real*8 R, beta, gamma, V, w, wt
        real*8 total, two_body, three_body, four_body
        real*8 r_k, r_l, r_m, x_He, y_He, z_He, xx, yy, zz
        real*8 v2, v3, v3_1, v3_2, v3_2_1, v3_2_2, v3_2_3, v3_2_4, v4, v4_1, v4_2, v4_3, v4_3_1, v4_3_2
        logical :: first = .false.

        if(.not.first) then
                call init
                print*, "Initialization of PES Done."
                first = .true.
        end if
		
        x_He = R*dsin(beta)*dcos(gamma)
        y_He = R*dsin(beta)*dsin(gamma)
        z_He = R*dcos(beta)

        two_body = 0.0d0
        do k = 1, 6
        xx = (x_He - x(k))**2d0
        yy = (y_He - y(k))**2d0
        zz = (z_He - z(k))**2d0
        r_k = dsqrt(dble(xx + yy + bz*zz))
        if(r_k.ne.r_k) then
                print*, R, beta, gamma
                stop
        end if

        v2 = w(r_k)**2d0 + c3*w(r_k)**3d0 + c4*w(r_k)**4d0 + c5*w(r_k)**5d0 + c6*wt(r_k)**6d0
        two_body = two_body + v2
        end do		

        three_body = 0.0d0
        do k = 1, 6
        xx = (x_He - x(k))**2d0
        yy = (y_He - y(k))**2d0
        zz = (z_He - z(k))**2d0
        r_k = dsqrt(dble(xx + yy + bz*zz))
        do l = 1, k-1
        xx = (x_He - x(l))**2d0
        yy = (y_He - y(l))**2d0
        zz = (z_He - z(l))**2d0
        r_l = dsqrt(dble(xx + yy + bz*zz))

        v3_1 = c11*w(r_k)*w(r_l) + c22*w(r_k)**2d0*w(r_l)**2d0
        v3_2_1 = c12*(w(r_k)*w(r_l)**2d0 + w(r_l)*w(r_k)**2d0)
        v3_2_2 = c13*(w(r_k)*w(r_l)**3d0 + w(r_l)*w(r_k)**3d0)
        v3_2_3 = c14*(w(r_k)*w(r_l)**4d0 + w(r_l)*w(r_k)**4d0)
        v3_2_4 = c23*(w(r_k)**2d0*w(r_l)**3d0 + w(r_l)**2d0*w(r_k)**3d0)
        v3_2 = v3_2_1 + v3_2_2 + v3_2_3 + v3_2_4
        v3 = v3_1 + v3_2
        three_body = three_body + v3
        end do
        end do

        four_body = 0.0d0
        do k = 1, 6
        xx = (x_He - x(k))**2d0
        yy = (y_He - y(k))**2d0
        zz = (z_He - z(k))**2d0
        r_k = dsqrt(dble(xx + yy + bz*zz))
        do l = 1, k-1
        xx = (x_He - x(l))**2d0
        yy = (y_He - y(l))**2d0
        zz = (z_He - z(l))**2d0
        r_l = dsqrt(dble(xx + yy + bz*zz))
        do m = 1, l-1
        xx = (x_He - x(m))**2d0
        yy = (y_He - y(m))**2d0
        zz = (z_He - z(m))**2d0
        r_m = dsqrt(dble(xx + yy + bz*zz))
        v4_1 = c111*w(r_k)*w(r_l)*w(r_m)
        v4_2 = c122*(w(r_k)*w(r_l)**2d0*w(r_m)**2d0 &
                + w(r_k)**2d0*w(r_l)*w(r_m)**2d0 &
                + w(r_k)**2d0*w(r_l)**2d0*w(r_m)) 
        v4_3_1 = c112*(w(r_k)*w(r_l)*w(r_m)**2d0 &
                + w(r_k)*w(r_l)**2d0*w(r_m) &
                + w(r_k)**2d0*w(r_l)*w(r_m))
        v4_3_2 = c113*(w(r_k)*w(r_l)*w(r_m)**3d0 &
                + w(r_k)*w(r_l)**3d0*w(r_m) &
                + w(r_k)**3d0*w(r_l)*w(r_m))
        v4_3 = v4_3_1 + v4_3_2
        v4 = v4_1 + v4_2 + v4_3
        four_body = four_body + v4
        end do
        end do
        end do

        total = two_body + three_body + four_body
        V = c0 + W0*total - shift_bikram
        return
end subroutine

function w(rr)
        use params
        implicit none
        real*8 rr, w
        w = 1.0d0 - exp(-a*(rr - r0))
        return
end function

function wt(rr)
        use params
        implicit none
        real*8 rr, wt
        if(rr.ge.r0) then
                wt = 1.0d0 - exp(-a*(rr - r0))
        else if(rr.lt.r0) then
                wt = 0.0d0
        else
                print*, "Something is wrong in W_T, Please Check", rr
                stop
        end if
        return
end function


