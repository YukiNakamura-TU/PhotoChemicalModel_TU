program main

  implicit none

  integer, parameter :: ni = 96
  integer :: i, ialt
  real(8) :: n_CO2(1:ni), n_N2(1:ni), n_O(1:ni), n_H(1:ni)
  real(8) :: t_n(1:ni), t_i(1:ni), t_e(1:ni)
  real(8) :: p_CO2i(1:ni), p_O2i(1:ni), p_Oi(1:ni), p_Hi(1:ni)
  real(8) :: altitude(1:ni)
  real(8) :: a1, a2, a3, a4, a5
  character :: c

  open(11, file = "min_neutral.dat", form ="formatted", status = "old")
  do i = 1, 6
    read(11, *) c
  end do
  do i = 1, ni
    read(11, *)  ialt, a1, a2, a3, a4, a5
    altitude(i) = real(ialt)
    n_CO2(i) = 1.0e+6 * a1
    n_N2(i) = 1.0e+6 * a2
    n_O(i) = 1.0e+6 * (a3 + a5)
    n_H(i) = 1.0e+6 * a4
  end do
  close(11)

  open(21, file = "CO2.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), n_CO2(i)
  end do
  close(21)
  open(21, file = "N2.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), n_N2(i)
  end do
  close(21)
  open(21, file = "O.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), n_O(i)
  end do
  close(21)
  open(21, file = "H.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), n_H(i)
  end do
  close(21)

  open(11, file = "min_temperature.dat", form ="formatted", status = "old")
  do i = 1, 6
    read(11, *) c
  end do
  do i = 1, ni
    read(11, *)  ialt, a1, a2, a3
    altitude(i) = real(ialt)
    t_n(i) = a1
    t_i(i) = a2
    t_e(i) = a3
  end do
  close(11)
  open(21, file = "T_e.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), t_e(i)
  end do
  close(21)
  open(21, file = "T_i.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), t_i(i)
  end do
  close(21)
  open(21, file = "T_n.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), t_n(i)
  end do
  close(21)

  open(11, file = "min_photoi.dat", form ="formatted", status = "old")
  do i = 1, 6
    read(11, *) c
  end do
  do i = 1, ni
    read(11, *)  ialt, a1, a2, a3,a4
    altitude(i) = real(ialt)
    p_CO2i(i) = 1.0e+6 * a1
    p_O2i(i) = 1.0e+6 * a2
    p_Oi(i) = 1.0e+6 * a3
    p_Hi(i) = 1.0e+6 * a4
  end do
  close(11)
  open(21, file = "photoi.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), p_CO2i(i), p_O2i(i), p_Oi(i), p_Hi(i)
  end do
  close(21)

end program main
