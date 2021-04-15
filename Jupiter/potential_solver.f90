program jupiter
 implicit none

!parameter

    double precision, parameter :: pi = dacos(-1.0d0)    ! pi
    double precision :: mu0 = 4.0d0 * pi * 1.0d-7 !magnetic permeability [H/m]

    integer i    ! a grid number of longitude directing eastward.
                 ! i = 1 is at midnight.
    integer j    ! a grid number of latitude directing northward.
                 ! j = 1 is at the south pole.
    integer h   !Tao's model altitude number.

    integer, parameter :: tn_i = 361    !the total number of i
    integer, parameter :: tn_j = 181    !the total number of j
    integer, parameter :: tn_h = 61    !the total number of h

    integer s   !iteration number
    integer n   !Hill current when n=1, R2 + Hill when n=2

    !the grid interval [rad]
    double precision :: delta_p = 2.0d0 * pi / (dble(tn_i) - 1.0d0)
    double precision :: delta_t = pi / (dble(tn_j) - 1.0d0)

    !(a radius of a planet + 2000[km]) * -1.0 [m]
    double precision :: R0 = - 7.3492d+7

    !magnetic moment [T m**3]
    double precision :: m(3) = (/0.0d0, 0.0d0, 1.5639d+20/)

    double precision longitude(tn_i)  !directing eastward. 0~360[degree]
    double precision latitude(tn_j)   !-90~90[degree]
    double precision phi_long(tn_i)   !longitude [rad]
    double precision theta(tn_j)      !latitude at (i,j) [rad]

    !variables of dip angle
    double precision  r(tn_i,tn_j,3)     ! a position vector from the planet center

    !dip angle I
    double precision cosI(tn_i,tn_j)    !cosine
    double precision sinI(tn_i,tn_j)    !sine
    double precision sinI_tmp(tn_i,tn_j,tn_h), cosI_tmp(tn_i,tn_j,tn_h)

    !variables of conductivity
    double precision hic_tt(tn_i,tn_j) !theta-theta
    double precision hic_pp(tn_i,tn_j) !phi-phi
    double precision hic_tp(tn_i,tn_j) !theta-phi

    double precision hic_tt_t(tn_i,tn_j) !theta  differential
    double precision hic_tt_p(tn_i,tn_j) !phi differential
    double precision hic_pp_t(tn_i,tn_j) !theta differential
    double precision hic_pp_p(tn_i,tn_j) !phi differential
    double precision hic_tp_t(tn_i,tn_j) !theta differential
    double precision hic_tp_p(tn_i,tn_j) !phi differential

    !field aligned current [A/m2]. downward is positive.
    double precision TotalFAC(tn_i,tn_j)      !Global input electric current
    double precision R2FAC(tn_i,tn_j)    !Region 2 current
    double precision HillC(tn_i,tn_j)    !Hill current
    double precision FAC(tn_i,tn_j)

    !Ionospheric electric potential
    double precision phi_new(tn_i,tn_j)
    double precision phi_old(tn_i,tn_j)

    double precision phi_dif, phi_sum, phi_tmp1, phi_tmp2, tmp(tn_i,tn_j)

    !Magnetospheric potential
    double precision phi_Hill(tn_i,tn_j) !Hill current only
    double precision phi_R2(tn_i,tn_j)   !Hill + R2
    double precision phi_cr(tn_i, tn_j)  !corotation
    double precision phi_scr(tn_i,tn_j)  !sub-corotation
    double precision phi_act(tn_i,tn_j)  !actual
    double precision phi_dd(tn_i,tn_j)   !by dawn-to-dusk E

    !dusk side Region2 corrent intensity (+:inward flow, -:outward flow)
        !North
    double precision :: r2N_long = 90.0d0 !longitude of peak [degree]
    double precision :: r2N_lat = 77.0d0  !latitude of peak [degree]
    double precision    r2N_intensity        !peak intensity
    double precision :: deltaN_lat = 4.0d0   !e-folding distance of lat. [degree]
    double precision :: deltaN_long = 45.0d0 !e-folding distance of long. [degree]

    double precision denominator1, denominator2
    double precision dlong1plus(tn_i,tn_j), dlong1minus(tn_i,tn_j)
    double precision dlong2plus(tn_i,tn_j), dlong2minus(tn_i,tn_j)

    !variables of a bird's eye view
    double precision rv(tn_j)    ! a radius vector from the polars

    !L shell (L value)
    double precision Lv(tn_j)

    !Electrical Field on equatorial plane
        !phi_long: directing longitude, Lv:directing L value
    double precision E_phi_long(tn_i,tn_j)
    double precision E_Lv(tn_i,tn_j)
        !x: directing anti-Sun ward, y: directing the orbital motion
    double precision E_y_Hill(tn_i,tn_j)
    double precision E_x_Hill(tn_i,tn_j)
    double precision E_y_R2(tn_i,tn_j)
    double precision E_x_R2(tn_i,tn_j)
    double precision E_y_dd(tn_i,tn_j)
    double precision E_x_dd(tn_i,tn_j)

    !Electric field on dawn-dusk cross section
    double precision E_Hill_dawn(tn_j)
    double precision E_Hill_dusk(tn_j)
    double precision E_R2_dawn(tn_j)
    double precision E_R2_dusk(tn_j)
    double precision E_cr_dawn(tn_j)
    double precision E_cr_dusk(tn_j)
    double precision E_scr_dawn(tn_j)
    double precision E_scr_dusk(tn_j)
    double precision E_act_dawn(tn_j)
    double precision E_act_dusk(tn_j)
    double precision E_dd_dawn(tn_j)
    double precision E_dd_dusk(tn_j)

    !sub-corotational potential at 5.9RJ
    double precision phi_cr_IPT_dawn
    double precision phi_cr_IPT_dusk

    !L shell of IPT after shift
    double precision Lv_shift_scr_dawn
    double precision Lv_shift_scr_dusk
    double precision Lv_shift_act_dawn
    double precision Lv_shift_act_dusk

    !System III rotaiton period [s] (= 9h55m29.37s)
    double precision :: T_rot = 35729.37d0

    !others
    double precision tmp1, tmp2, tmp3, tmp4
    double precision tmp5, tmp6, tmp7, tmp8

!reading
    open(11, file = './FAC_R2_60d6.dat', status = 'unknown')
        do j = 1, tn_j
        do i = 1, tn_i
            read(11,*) R2FAC(i,j), HillC(i,j), TotalFAC(i,j)
        end do
        end do
    close(11)

    open(12, file = './dipole.dat', status = 'old')
        do h = 1, tn_h
        do j = 1, tn_j
        do i = 1, tn_i
            read(12,*) tmp1, sinI_tmp(i,j,h), cosI_tmp(i,j,h)
        end do
        end do
        end do
    close(12)

    !We use I at altitude 1000 [km] as height-integrated calculation
    do j = 1, tn_j
    do i = 1, tn_i
        sinI(i,j) = sinI_tmp(i,j,20)
        cosI(i,j) = cosI_tmp(i,j,20)
    end do
    end do

!coordinate
    do i = 1, tn_i
        longitude(i) = dble(i) - 1.0d0
        phi_long(i) = longitude(i) * pi / 180.0d0
    end do

    do j = 1, tn_j
        latitude(j) = dble(j) - 91.0d0
        theta(j) = latitude(j) * pi / 180.0d0
    end do

    do j = 1, tn_j
    do i = 1, tn_i
        r(i,j,1) = dabs(R0) * dcos(theta(j)) * dcos(phi_long(i))
        r(i,j,2) = dabs(R0) * dcos(theta(j)) * dsin(phi_long(i))
        r(i,j,3) = dabs(R0) * dsin(theta(j))
    end do
    end do

!Two cases of "only Hill current" and "Hill + R2 current"
do n = 1, 2

    if( n == 1 ) then
        open(13, file = './parallel/cond_Metal_60MA/hic_Hill.dat', status = 'unknown' )
            do j = 1, tn_j
            do i = 1, tn_i
                read(13,*) hic_tt(i,j), hic_pp(i,j), hic_tp(i,j)
            end do
            end do
        close(13)
    else if( n == 2) then
        open(14, file = './parallel/cond_Metal_60MA/hic_H_R2.dat', status = 'unknown' )
            do j = 1, tn_j
            do i = 1, tn_i
                read(14,*) hic_tt(i,j), hic_pp(i,j), hic_tp(i,j)
            end do
            end do
        close(14)
    end if

    if( n == 1 ) then
        FAC(:,:) = HillC(:,:)   !(:,:) = All (i,j)
    else if( n == 2 ) then
        FAC(:,:) = TotalFAC(:,:)
    end if

    !a gradient of hic_pp, hic_tt, hic_tp
    do j = 2, 180
    do i = 2, 360
        hic_pp_p(i,j) = ( hic_pp(i+1,j) - hic_pp(i-1,j) ) * 0.5d0 / delta_p
        hic_pp_t(i,j) = ( hic_pp(i,j+1) - hic_pp(i,j-1) ) * 0.5d0 / delta_t
        hic_tt_p(i,j) = ( hic_tt(i+1,j) - hic_tt(i-1,j) ) * 0.5d0 / delta_p
        hic_tt_t(i,j) = ( hic_tt(i,j+1) - hic_tt(i,j-1) ) * 0.5d0 / delta_t
        hic_tp_p(i,j) = ( hic_tp(i+1,j) - hic_tp(i-1,j) ) * 0.5d0 / delta_p
        hic_tp_t(i,j) = ( hic_tp(i,j+1) - hic_tp(i,j-1) ) * 0.5d0 / delta_t
    end do
    end do

    do j = 2, 180
        hic_pp_p(1,j) = ( hic_pp(2,j) - hic_pp(360,j) ) * 0.5d0 / delta_p
        hic_pp_p(361,j) = hic_pp_p(1,j)
        hic_pp_t(1,j) = ( hic_pp(1,j+1) - hic_pp(1,j-1) ) * 0.5d0 / delta_t
        hic_pp_t(361,j) = hic_pp_t(1,j)
        hic_tt_p(1,j) = ( hic_tt(2,j) - hic_tt(360,j) ) * 0.5d0 / delta_p
        hic_tt_p(361,j) = hic_tt_p(1,j)
        hic_tt_t(1,j) = ( hic_tt(1,j+1) - hic_tt(1,j-1) ) * 0.5d0 / delta_t
        hic_tt_t(361,j) = hic_tt_t(1,j)
        hic_tp_p(1,j) = ( hic_tp(2,j) - hic_tp(360,j) ) * 0.5d0 / delta_p
        hic_tp_p(361,j) = hic_tp_p(1,j)
        hic_tp_t(1,j) = ( hic_tp(1,j+1) - hic_tp(1,j-1) ) * 0.5d0 / delta_t
        hic_tp_t(361,j) = hic_tp_t(1,j)
    end do

    do i = 2, 360
        hic_pp_p(i,1) = 0.0d0
        hic_pp_p(i,181) = 0.0d0
        hic_pp_t(i,1) = 0.0d0
        hic_pp_t(i,181) = 0.0d0
        hic_tt_p(i,1) = 0.0d0
        hic_tt_p(i,181) = 0.0d0
        hic_tt_t(i,1) = 0.0d0
        hic_tt_t(i,181) = 0.0d0
        hic_tp_p(i,1) = 0.0d0
        hic_tp_p(i,181) = 0.0d0
        hic_tp_t(i,1) = 0.0d0
        hic_tp_t(i,181) = 0.0d0
    end do

    !ionospheric potential
    phi_old(1:tn_i, 1:tn_j) = 0.0d0

    do s = 1, 100000
        do i = 2, 360
        do j = 2, 180

            phi_new(i,j) &
&           = (R0 * cos(theta(j)) * delta_t * delta_p)**2.0d0 * 0.5d0 &
&           / (cos(theta(j))**2.0d0 * hic_tt(i,j) * delta_p**2.0d0 + hic_pp(i,j) &
&           * delta_t**2.0d0) * ( FAC(i,j) * abs(sinI(i,j)) &
&           + ((-sin(theta(j)) / cos(theta(j)) * hic_tt(i,j) + hic_tt_t(i,j) &
&           - hic_tp_p(i,j) / cos(theta(j))) * 0.5d0 / R0**2.0d0 / delta_t &
&           + hic_tt(i,j) / R0**2.0d0 / delta_t**2.0d0) * phi_old(i,j+1) &
&           + ((sin(theta(j)) / cos(theta(j)) * hic_tt(i,j) - hic_tt_t(i,j) &
&           + hic_tp_p(i,j) / cos(theta(j)) ) * 0.5d0 / R0**2.0d0 / delta_t &
&           + hic_tt(i,j) / R0**2.0d0 / delta_t**2.0d0) * phi_old(i,j-1) &
&           + ((hic_tp_t(i,j) / cos(theta(j)) + hic_pp_p(i,j) / cos(theta(j))**2.0d0) &
&           * 0.5d0 / R0**2.0d0 / delta_p + hic_pp(i,j) / R0**2.0d0 &
&           / cos(theta(j))**2.0d0 / delta_p**2.0d0) * phi_old(i+1,j) &
&           + ((-hic_tp_t(i,j) / cos(theta(j)) - hic_pp_p(i,j) / cos(theta(j))**2.0d0)&
&           * 0.5d0 / R0**2.0d0 / delta_p + hic_pp(i,j) / R0**2.0d0 &
&           / cos(theta(j))**2.0d0 / delta_p**2.0d0) * phi_old(i-1,j) )

            phi_new(1,j) &
&           = (R0 * cos(theta(j)) * delta_t * delta_p)**2.0d0 * 0.5d0 &
&           / (cos(theta(j))**2.0d0 * hic_tt(1,j) * delta_p**2.0d0 + hic_pp(1,j) &
&           * delta_t**2.0d0) * ( FAC(1,j) * abs(sinI(1,j)) &
&           + ((-sin(theta(j)) / cos(theta(j)) * hic_tt(1,j) + hic_tt_t(1,j) &
&           - hic_tp_p(1,j) / cos(theta(j))) * 0.5d0 / R0**2.0d0 / delta_t &
&           + hic_tt(1,j) / R0**2.0d0 / delta_t**2.0d0) * phi_old(1,j+1) &
&           + ((sin(theta(j)) / cos(theta(j)) * hic_tt(1,j) - hic_tt_t(1,j) &
&           + hic_tp_p(1,j) / cos(theta(j))) * 0.5d0 / R0**2.0d0 / delta_t &
&           + hic_tt(1,j) / R0**2.0d0 / delta_t**2.0d0) * phi_old(1,j-1) &
&           + ((hic_tp_t(1,j) / cos(theta(j)) + hic_pp_p(1,j) / cos(theta(j))**2.0d0) &
&           * 0.5d0 / R0**2.0d0 / delta_p + hic_pp(1,j) / R0**2.0d0 &
&           / cos(theta(j))**2.0d0 / delta_p**2.0d0) * phi_old(2,j) &
&           + ((-hic_tp_t(1,j) / cos(theta(j)) - hic_pp_p(1,j) / cos(theta(j))**2.0d0)&
&           * 0.5d0 / R0**2.0d0 / delta_p + hic_pp(1,j) / R0**2.0d0 &
&           / cos(theta(j))**2.0d0 / delta_p**2.0d0) * phi_old(360,j) )

            phi_new(361,j) = phi_new(1,j)
        end do
        end do

        phi_tmp1 = sum( phi_new(1:360,2) ) / 360.0d0
        phi_tmp2 = sum( phi_new(1:360,180) ) / 360.0d0

        do i = 1, 361
            phi_new(i,1) = phi_tmp1
            phi_new(i,tn_j) = phi_tmp2
        end do

        do j = 1, tn_j
        do i = 1, tn_i
            if( isnan( phi_new(i,j) ) ) then
                open(200, file = './parallel/Progress_Phi.dat', status = 'unknown' )
                    write(200,*) 'error', i, j
                close(200)
            end if
        end do
        end do

        phi_dif = sum(dabs(phi_new(1:360, 1:tn_j) - phi_old(1:360, 1:tn_j)))
        phi_sum = sum(dabs(phi_new(1:360, 1:tn_j)))
        if ( mod(s, 1000) == 0 ) then
            open(200, file = './parallel/Progress_Phi.dat', status = 'unknown' )
                write(200, *) 'n =', n
                write(200, *) s, phi_dif, phi_sum, phi_dif / phi_sum
            close(200)
        end if
        if(phi_dif / phi_sum < 1.0d-5) exit

        phi_old(1:tn_i, 1:tn_j) = phi_new(1:tn_i, 1:tn_j)
    end do  !of s

    if( n == 1 ) then
            phi_Hill(:,:) = phi_new(:,:)    !(:,:) is (1:tn_i, 1:tn_j).
    else if ( n == 2 ) then
            phi_R2(:,:) = phi_new(:,:)
    end if

end do !of n

open(200, file = './parallel/Progress_Phi.dat', status = 'unknown' )
    write(200, *) 'done!'
close(200)

!L-shell
    do j = 1, tn_j-1
        Lv(j) = 1 / ( dcos(theta(j))**2.0d0 )
    end do

!Electrical potential by Hill current at equatorial plane
    do j = 92, tn_j-2
        do i = 2, 360
            E_phi_long(i,j) = -( phi_Hill(i+1,j) - phi_Hill(i-1,j) ) * 0.5d0 &
&                           / ( delta_p * Lv(j) * dabs(R0) )
            E_Lv(i,j) = -( phi_Hill(i,j+1) - phi_Hill(i,j-1) ) &
&                     / ( (Lv(j+1) - Lv(j-1)) *  dabs(R0) )
        end do
        E_phi_long(1,j) = -( phi_Hill(2,j) - phi_Hill(360,j) ) * 0.5d0 &
&                       / ( delta_p *  Lv(j) * dabs(R0) )
        E_phi_long(361,j) = E_phi_long(1,j)
    end do

    do j = 92, tn_j-2
    do i = 1, tn_i
        E_x_Hill(i,j) = E_Lv(i,j) * dcos(phi_long(i)) &
&                     - E_phi_long(i,j) * dsin(phi_long(i))
        E_y_Hill(i,j) = E_Lv(i,j) * dsin(phi_long(i)) &
&                     + E_phi_long(i,j) * dcos(phi_long(i))
    end do
    end do

!Electrical potential by Hill + R2 at equatorial plane
    do j = 92, tn_j-2
        do i = 2, 360
            E_phi_long(i,j) = -( phi_R2(i+1,j) - phi_R2(i-1,j) ) * 0.5d0 &
&                           / ( delta_p * Lv(j) * dabs(R0) )
            E_Lv(i,j) = -( phi_R2(i,j+1) - phi_R2(i,j-1) ) &
&                     / ( (Lv(j+1) - Lv(j-1)) *  dabs(R0) )
        end do
        E_phi_long(1,j) = -( phi_R2(2,j) - phi_R2(360,j) ) * 0.5d0 &
&                       / ( delta_p *  Lv(j) * dabs(R0) )
        E_phi_long(361,j) = E_phi_long(1,j)
    end do

    do j = 92, tn_j-2
    do i = 1, tn_i
        E_x_R2(i,j) = E_Lv(i,j) * dcos(phi_long(i)) &
&                   - E_phi_long(i,j) * dsin(phi_long(i))
        E_y_R2(i,j) = E_Lv(i,j) * dsin(phi_long(i)) &
&                   + E_phi_long(i,j) * dcos(phi_long(i))
    end do
    end do

!Electrical field and potential by Hill current at dawn-dusk cross section
    do j = 92, tn_j-2
        E_Hill_dawn(j) = -( phi_Hill(91,j-1) - phi_Hill(91,j+1) ) &
&                      / ( abs(Lv(j+1) - Lv(j-1)) *  abs(R0) )
        E_Hill_dusk(j) = -( phi_Hill(271,j+1) - phi_Hill(271,j-1) ) &
&                      / ( abs(Lv(j+1) - Lv(j-1)) *  abs(R0) )
    end do

!Electric field and potential by Hill and R2 current at dawn-dusk cross section
    do j = 91, tn_j-2
        E_R2_dawn(j) = -( phi_R2(91,j-1) - phi_R2(91,j+1) ) &
&                    / ( abs(Lv(j+1) - Lv(j-1)) *  abs(R0) )
        E_R2_dusk(j) = -( phi_R2(271,j+1) - phi_R2(271,j-1) ) &
&                    / ( abs(Lv(j+1) - Lv(j-1)) *  abs(R0) )
    end do

!Corotation Electric field and Potential
    do j = 91, 180
        E_cr_dawn(j) = -2.0d0 * pi * dabs(m(3)) &
&                    / ( T_rot * (Lv(j) * dabs(R0))**2.0d0 )
        E_cr_dusk(j) = 2.0d0 * pi * dabs(m(3)) &
&                    / ( T_rot * (Lv(j) * dabs(R0))**2.0d0 )
    end do

    do j = 91, 180
    do i = 1, tn_i
        phi_cr(i,j) = 2.0d0 * pi * dabs(m(3)) / ( T_rot * Lv(j) * dabs(R0) )
    end do
    end do

!Sub-Corotetion Electric field and potential
    do j = 91, tn_j-1
    do i = 1, tn_i
        phi_scr(i,j) = phi_cr(i,j) + phi_Hill(i,j)
    end do
    end do

    do j = 91, tn_j-2
        E_scr_dawn(j) = -( phi_scr(91,j-1) - phi_scr(91,j+1) ) &
&                    / ( dabs(Lv(j+1) - Lv(j-1)) *  dabs(R0) )
        E_scr_dusk(j) = -( phi_scr(271,j+1) - phi_scr(271,j-1) ) &
&                    / ( dabs(Lv(j+1) - Lv(j-1)) *  dabs(R0) )
    end do

!Actual electric field and potential
    do j = 91, tn_j-1
    do i = 1, tn_i
        phi_act(i,j) = phi_cr(i,j) + phi_R2(i,j)
    end do
    end do

    do j = 91, tn_j-2
        E_act_dawn(j) = -( phi_act(91,j-1) - phi_act(91,j+1) ) &
&                    / ( dabs(Lv(j+1) - Lv(j-1)) *  dabs(R0) )
        E_act_dusk(j) = -( phi_act(271,j+1) - phi_act(271,j-1) ) &
&                    / ( dabs(Lv(j+1) - Lv(j-1)) *  dabs(R0) )
    end do

!Dawn-to-Dusk electric field and potential
    do j = 91, tn_j-1
    do i = 1, tn_i
        phi_dd(i,j) = phi_act(i,j) - phi_scr(i,j)
    end do
    end do

    do j = 91, tn_j-2
        E_dd_dawn(j) = -( phi_dd(91,j-1) - phi_dd(91,j+1) ) &
&                    / ( dabs(Lv(j+1) - Lv(j-1)) *  dabs(R0) )
        E_dd_dusk(j) = -( phi_dd(271,j+1) - phi_dd(271,j-1) ) &
&                    / ( dabs(Lv(j+1) - Lv(j-1)) *  dabs(R0) )
    end do

!Dawn-to-Dusk Electrical potential of LT variation
    do j = 92, tn_j-2
        do i = 2, 360
            E_phi_long(i,j) = -( phi_dd(i+1,j) - phi_dd(i-1,j) ) * 0.5d0 &
&                           / ( delta_p * Lv(j) * dabs(R0) )
            E_Lv(i,j) = -( phi_dd(i,j+1) - phi_dd(i,j-1) ) &
&                     / ( dabs(Lv(j+1) - Lv(j-1)) *  dabs(R0) )
        end do
        E_phi_long(1,j) = -( phi_dd(2,j) - phi_dd(360,j) ) * 0.5d0 &
&                       / ( delta_p *  Lv(j) * dabs(R0) )
        E_phi_long(361,j) = E_phi_long(1,j)
    end do

    do j = 92, tn_j-2
    do i = 1, tn_i
        E_x_dd(i,j) = E_Lv(i,j) * dcos(phi_long(i)) &
&                   - E_phi_long(i,j) * dsin(phi_long(i))
        E_y_dd(i,j) = E_Lv(i,j) * dsin(phi_long(i)) &
&                   + E_phi_long(i,j) * dcos(phi_long(i))
    end do
    end do

!IPT shift amount at dawn side
    do j = 91, tn_j-1
        if( Lv(j) <= 5.9d0 .and. 5.9d0 < Lv(j+1) ) then
            phi_cr_IPT_dawn &
&           = phi_cr(91,j) + ( phi_cr(91,j+1) - phi_cr(91,j) ) &
&           / ( Lv(j+1) - Lv(j) ) * ( 5.9d0 - Lv(j) )

            exit
        end if
    end do

    do j = 91, tn_j-1
        if( phi_scr(91,j+1) <= phi_cr_IPT_dawn &
&           .and. phi_cr_IPT_dawn < phi_scr(91,j) ) then
            Lv_shift_scr_dawn &
&           = Lv(j) + ( Lv(j+1) - Lv(j) ) / ( phi_scr(91,j+1) - phi_scr(91,j) ) &
&           * ( phi_cr_IPT_dawn - phi_scr(91,j) )

            exit
        end if
    end do

    do j = 91, tn_j-1
        if( phi_act(91,j+1) <= phi_cr_IPT_dawn &
&           .and. phi_cr_IPT_dawn < phi_act(91,j) ) then
            Lv_shift_act_dawn &
&           = Lv(j) + ( Lv(j+1) - Lv(j) ) / ( phi_act(91,j+1) - phi_act(91,j) ) &
&           * ( phi_cr_IPT_dawn - phi_act(91,j) )

            exit
        end if
    end do

!IPT shift amount at dusk side
    do j = 91, tn_j-1
        if( Lv(j) <= 5.9d0 .and. 5.9d0 < Lv(j+1) ) then
            phi_cr_IPT_dusk &
&           = phi_cr(271,j) + ( phi_cr(271,j+1) - phi_cr(271,j) ) &
&           / ( Lv(j+1) - Lv(j) ) * ( 5.9d0 - Lv(j) )

            exit
        end if
    end do

    do j = 91, tn_j-1
        if( phi_scr(271,j+1) <= phi_cr_IPT_dusk &
&           .and. phi_cr_IPT_dusk < phi_scr(271,j) ) then
            Lv_shift_scr_dusk &
&           = Lv(j) + ( Lv(j+1) - Lv(j) ) / ( phi_scr(271,j+1) - phi_scr(271,j) ) &
&           * ( phi_cr_IPT_dusk - phi_scr(271,j) )

            exit
        end if
    end do

    do j = 91, tn_j-1
        if( phi_act(271,j+1) <= phi_cr_IPT_dusk &
&           .and. phi_cr_IPT_dusk < phi_act(271,j) ) then
            Lv_shift_act_dusk &
&           = Lv(j) + ( Lv(j+1) - Lv(j) ) / ( phi_act(271,j+1) - phi_act(271,j) ) &
&           * ( phi_cr_IPT_dusk - phi_act(271,j) )

            exit
        end if
    end do

!bird's eye view
    do j = 1, tn_j
        rv(j) = dabs(R0) * ( 0.5d0 * pi - dabs(theta(j)) )
    end do

!print
    open(20, file = './parallel/Potential_Metal_60MA/Ionospheric_potential_map.dat', status = 'unknown')
        do j = 1, 91
        do i = 1, tn_i
            write(20, *) phi_long(i), phi_Hill(i,j+90)*1.0d-3, &
&                        phi_R2(i,j+90)*1.0d-3, rv(j+90), &
                         phi_Hill(i,j)*1.0d-3, phi_R2(i,j)*1.0d-3, rv(j)
        end do
            write(20, *)
        end do
    close(20)

    open(22, file = './parallel/Potential_Metal_60MA/phi_dawndusk.dat', status = 'unknown')
        do j = tn_j-1, 91, -1
            if( Lv(j) <= 10 ) then
                write(22,*) -Lv(j), phi_cr(91,j)*1.0d-3, &
&                           phi_Hill(91,j)*1.0d-3, phi_scr(91,j)*1.0d-3, &
&                           phi_R2(91,j)*1.0d-3, phi_dd(91,j)*1.0d-3, &
&                           phi_act(91,j)*1.0d-3
            end if
        end do
        do j = 91, tn_j-1
            if( Lv(j) <= 10 ) then
                write(22,*) Lv(j), phi_cr(271,j)*1.0d-3, &
&                           phi_Hill(271,j)*1.0d-3, phi_scr(271,j)*1.0d-3, &
&                           phi_R2(271,j)*1.0d-3, phi_dd(271,j)*1.0d-3, &
&                           phi_act(271,j)*1.0d-3
            end if
        end do
    close(22)

    open(23, file = './parallel/Potential_Metal_60MA/E_dawndusk.dat', status = 'unknown')
        do j = tn_j-1, 91, -1
            if( Lv(j) <= 10 ) then
                write(23,*) -Lv(j), E_cr_dawn(j)*1.0d3, &
&                           E_Hill_dawn(j)*1.0d3, E_scr_dawn(j)*1.0d3, &
&                           E_R2_dawn(j)*1.0d3, E_act_dawn(j)*1.0d3, &
&                           E_dd_dawn(j)*1.0d3
            end if
        end do
        do j = 91, tn_j-1
            if( Lv(j) <= 10 ) then
                write(23,*) Lv(j), E_cr_dusk(j)*1.0d3, &
&                           E_Hill_dusk(j)*1.0d3, E_scr_dusk(j)*1.0d3, &
&                           E_R2_dusk(j)*1.0d3, E_act_dusk(j)*1.0d3, &
&                           E_dd_dusk(j)*1.0d3
            end if
        end do
    close(23)

    open(24, file = './parallel/Potential_Metal_60MA/Magnetospheric_potential_map.dat', status = 'unknown')
        do j = 91, tn_j-1
        do i = 1, tn_i
            write(24,*) phi_long(i), phi_cr(i,j)*1.0d-3, &
&                       phi_Hill(i,j)*1.0d-3, phi_scr(i,j)*1.0d-3, &
&                       phi_R2(i,j)*1.0d-3, phi_act(i,j)*1.0d-3, &
&                       phi_dd(i,j)*1.0d-3, Lv(j)
        end do
            write(24,*)
        end do
    close(24)

    open(24, file = './parallel/Potential_Metal_60MA/E_LT.dat', status = 'unknown')
        do i = 1, tn_i-1
            write(24,*) i-1, -E_y_dd(i,136)*1.0d+3, -E_y_dd(i,151)*1.0d+3, &
&                       -E_y_dd(i,157)*1.0d+3
        end do
    close(24)

    open(25, file = './parallel/Potential_Metal_60MA/IPT shift amount.dat', status = 'unknown')
            write(25,*) 'dawn'
            write(25,*) 'corotarional potential at 5.9RJ is', phi_cr_IPT_dawn
            write(25,*) 'subcorotational L shell of IPT is', Lv_shift_scr_dawn
            write(25,*) 'actual L shell of IPT is', Lv_shift_act_dawn
            write(25,*) 'dusk'
            write(25,*) 'corotarional potential at 5.9RJ is', phi_cr_IPT_dusk
            write(25,*) 'subcorotational L shell of IPT is', Lv_shift_scr_dusk
            write(25,*) 'actual L shell of IPT is', Lv_shift_act_dusk
    close(25)


end program jupiter
