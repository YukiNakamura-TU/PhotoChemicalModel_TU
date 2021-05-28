program main

  implicit none

  integer, parameter :: ni = 96
  integer, parameter :: ni_0 = 251
  integer :: i, i0
  real(8) :: n_O2(1:ni_0)
  real(8) :: altitude(1:ni_0)
  real(8) :: a1, a2

  i0 = 51
  open(11, file = "O2_nakamura.dat", form ="formatted", status = "old")
  do i = 1, ni_0
    read(11, *)  a1, a2
    altitude(i) = a1
    n_O2(i) = a2
  end do
  close(11)

  open(21, file = "./density/O2.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i0), n_O2(i0)
    i0 = i0 + 2
  end do
  close(21)

end program main
