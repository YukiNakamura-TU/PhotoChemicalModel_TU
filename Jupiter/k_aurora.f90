program Jupiter
implicit none

!parameter

    !phisical constant
    double precision :: pi = dacos(-1.0d0)

    integer ix, iy, iz, iv, i

    integer, parameter :: nx = 361    !the total number of i
    integer, parameter :: ny = 181    !the total number of j
    integer, parameter :: nz = 141     !the total number of h

    double precision theta(ny)        !latitude at (i,j) [rad]
    double precision phi_long(nx)     !longitude at (i,j) [rad]

    double precision alt(0:nz+1)  !altitude at (h) [m]
    double precision delta_alt(0:nz+1)  !the gap of h [m]

    !mass [kg]
    double precision :: m_e = 9.109d-31  !electron
    double precision :: m_H2 =  0.3347d-26

    !quantum of electricy [C]
    double precision :: q_e = 1.60217662d-19

    !Aurora electron Precipitations

    integer, parameter :: nv = 18751
    double precision v(nv) ![m/s]
    double precision :: delta_v = 1.0d5

    double precision v_0
    double precision :: W_th = 4.005d-16 ![J] = 2.5d3 [eV]
    double precision :: j_para_i_0 = 0.0134d-6 ![A/m^2]
    double precision :: j_para_i ![A/m^2]   !FAC
    double precision f_0
    double precision f(nv)
 
    double precision epsilon_ele(nv) !electron energy [eV]

    double precision flux(nv) ![m^-2 s^-1 eV^-1]

    !Ionization rate by aurora electron
    double precision k_edv(nv)    !the energy dependent variable part of lambda
    double precision lambda_0(nz,nv)
                            !the height dependent variable part of lambda
    double precision R_cmd(nz) !the column mass density above a corresponding alt.
    double precision R_0_cmd(nv)  !the column mass density
    double precision x_cmd(nz,nv) != R_cmd / R_0_cmd
    double precision tmp_ion(nz,nv)
    double precision eta_tilt(nz,nv)
    double precision a_tilt
    double precision b_tilt
    double precision c_tilt
    double precision q_ion(nz,nv)

    !the angle the magnetic field and the local vertical at h = nz [rad]
    double precision theta_tilt(nx,ny)

    !number density of fixed neutral spacies [m^-3]
    double precision n_H2(nz+1)  

    !column number density above alt(h)
    double precision cnd_H2(nz)

    !Rates [m^3 s^-1] (k1-k6 is production rate [m^-3 s^-1])
    double precision k1(nz)
    double precision k1arr(1000,nz)
    !others
    double precision tmp, tmp1, tmp2
    character filename*128

!latitude and longitude [rad]
    do ix = 1, nx
        phi_long(ix) = ( dble(ix) - 1.0d0 ) * pi / 180.0d0
    end do

    do iy = 1, ny
        theta(iy) = ( dble(iy) - 91.0d0 ) * pi / 180.0d0
    end do

    do iz = 0, nz+1
        delta_alt(iz) = 20.0d3
        alt(iz) = 200.0d3 + dble(iz-1) * delta_alt(iz)
    end do

!read

    open(14, file = './dipole_tilt.dat', status = 'old')
        do iy = 1, ny
        do ix = 1, nx
            read(14,*) theta_tilt(ix,iy)
        end do
        end do
    close(14)

    open(14, file = './metal_Hill/input/density/H2.dat', status = 'old')
        do iz = 1, nz
            read(14,*) tmp, n_H2(iz)
        end do
    close(14)

    cnd_H2 = 0.0d0
    do iz = nz, 1, -1
        cnd_H2(iz) = cnd_H2(iz) + n_H2(iz)*delta_alt(iz)
    end do


! calculating ionization rate by Tao et al., 2009
        ix = 181
        iy = 164 

        do i = 1, 100
            
            j_para_i = 1.0d-6 * dble(2*i)/100.0d0 ! 0.02-2.0 uA/m^2

            print *, 'i = ',i, '; j_para_i = ', j_para_i*1.0d6, 'uA/m^2'
        
                !Aourora electron precipitaion (Aurora's electron precipitates downward.)

                v_0 = dsqrt( 2.0d0 * W_th / m_e * (j_para_i / j_para_i_0 - 1.0d0) )

                do iv = 1, nv
                    v(iv) = 5.0d5 + dble(iv-1) * delta_v
                end do

                f_0 = j_para_i / ( pi * q_e * dsqrt(3.0d0)*pi/9.0d0 * v_0**4.0d0 )

                do iv = 1, nv
                    tmp = v(iv)/v_0
                    f(iv) = f_0 / ( tmp**2.0d0 + tmp**8.0d0 )
                end do

                do iv = 1, nv
                    epsilon_ele(iv) = 0.5d0 * m_e * v(iv)**2.0d0 / 1.602d-19
                end do

                do iv = 1, nv
                    flux(iv) = v(iv)**2.0d0 / m_e * f(iv) * 1.602d-19
                end do

                !Ionization rate by aurora precipitation
                do iv = 1, nv
                    k_edv(iv) = 0.93d0 - 0.444d0 &
            &                * dtanh( 0.8d0 * dlog10(epsilon_ele(iv)) - 3.38d0 )
                end do

                do iv = 1, nv
                    R_0_cmd(iv) = 3.31d-5 * ( epsilon_ele(iv) * 1.0d-3)**1.41d0
                end do

                do iz = 1, nz
                    R_cmd(iz) = m_H2 * cnd_H2(iz)
                end do

                do iv = 1, nv
                    do iz = 1, nz
                        x_cmd(iz,iv) = R_cmd(iz) / R_0_cmd(iv)
                    end do
                end do

                do iv = 1, nv
                    do iz = 1, nz
                        if ( 0.0d0 <= x_cmd(iz,iv) .and. x_cmd(iz,iv) <= 0.04d0 ) then
                            lambda_0(iz,iv) &
                &           = x_cmd(iz,iv) * (-578.9d0 * x_cmd(iz,iv) + 45.16d0) + 1.271d0
                        else if ( 0.04d0 < x_cmd(iz,iv) .and. x_cmd(iz,iv) <= 0.6d0 ) then
                            lambda_0(iz,iv) &
                &           = x_cmd(iz,iv) * ( x_cmd(iz,iv) * ( x_cmd(iz,iv) * ( x_cmd(iz,iv) &
                &           * ( x_cmd(iz,iv) * (-780.2d0 * x_cmd(iz,iv) + 1876.0d0) - 1819.0d0 ) &
                &           + 904.4d0 ) - 234.5d0 ) + 23.04d0 ) + 1.613d0
                        else if ( 0.6d0 < x_cmd(iz,iv) .and. x_cmd(iz,iv) <= 1.0d0 ) then
                            lambda_0(iz,iv) &
                &           = x_cmd(iz,iv) * ( x_cmd(iz,iv) * (-5.868d0 * x_cmd(iz,iv) + 15.12d0) &
                &           - 13.12d0 ) + 3.851d0
                        else if ( 1.0d0 < x_cmd(iz,iv) ) then
                            lambda_0(iz,iv) = 0.0d0
                        end if
                    end do
                end do

                c_tilt = - dlog10( theta_tilt(ix,iy) &
            &          * (-0.0436d0 * theta_tilt(ix,iy) + 0.0499) + 0.0302d0 )
                b_tilt = - dlog10( 0.75d0 * dcos(theta_tilt(ix,iy)) + 0.25d0 ) &
            &          / ( c_tilt - 0.523d0 )
                a_tilt = 2.0d0 * ( 1.0d0 - dcos(theta_tilt(ix,iy)) ) &
            &          / ( dcos(theta_tilt(ix,iy)) * (6.0d0 - c_tilt)**2.0d0 )

                tmp1 = 10.0d0**(-c_tilt)
                do iv = 1, nv
                    do iz = 1, nz
                        if( 0 <= x_cmd(iz,iv) .and. x_cmd(iz,iv) <= 1.0d-6 ) then
                            eta_tilt(iz,iv) = 1.0d0
                        else if( 1.0d-6 < x_cmd(iz,iv) .and. x_cmd(iz,iv) <= tmp1 ) then
                            tmp2 = dlog10( x_cmd(iz,iv) )
                            eta_tilt(iz,iv) &
                &           = -a_tilt * ( tmp2 + 6.0d0 ) * ( tmp2 + c_tilt ) + 1.0d0
                        else if( tmp1 < x_cmd(iz,iv) .and. x_cmd(iz,iv) <= 1.0d0 ) then
                            eta_tilt(iz,iv) = ( tmp1 / x_cmd(iz,iv) )**b_tilt
                        else if( 1.0d0 < x_cmd(iz,iv) ) then
                            eta_tilt(iz,iv) = 0.0d0
                        end if
                    end do
                end do

                do iv = 1, nv
                    do iz = 1, nz
                        q_ion(iz,iv) = epsilon_ele(iv) * k_edv(iv) * lambda_0(iz,iv) &
                &                  / ( 30.0d0 * R_0_cmd(iv) ) * m_H2 * n_H2(iz)

                        tmp_ion(iz,iv) = q_ion(iz,iv)* eta_tilt(iz,iv) * flux(iv) &
                &                    * m_e * v(iv) / 1.602d-19
                    end do
                end do

                do iz = 1, nz
                    tmp = 0.0d0
                    do iv = 1, nv-1
                        tmp = tmp + 0.5d0 * ( tmp_ion(iz,iv+1) + tmp_ion(iz,iv) ) * delta_v
                    end do
                    k1(iz) = tmp
                end do

                if ( j_para_i_0 <= j_para_i .and. j_para_i < 0.02d-6) then
                    do iz = 1, nz
                        k1(iz) = k1(iz) * 1.0d-1 + ( 0.9d0 * k1(iz) ) &
            &                     / ( 0.02d-6 - j_para_i_0 ) * ( j_para_i - j_para_i_0 )
                    end do
                else if ( 0.0d0 < j_para_i .and. j_para_i < j_para_i_0 ) then
                    do iz = 91, nz !!!!!!!!!!! what is this ?????????????????????????
                        k1(iz) = k1(iz) * 1.0d-1 !!!!!!!!!!! what is this ?????????????????????????
                    end do
                end if

                !write(*,*) 'j_para_i = ',j_para_i
                do iz = 1, nz
                    !if ( isnan(k1(iz)) ) then
                        !write(*,*) iz, k1(iz)/n_H2(iz)
                    !end if
                    k1arr(i,iz) = k1(iz)
                end do


        end do

        open(11, file = './k1.dat', status = 'unknown' )
            do iz = 1, nz
                write(11, *) alt(iz), (k1arr(i,iz), i = 1, 100)
            end do
        close(11)





!print



end program Jupiter
