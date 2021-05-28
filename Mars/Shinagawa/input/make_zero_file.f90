program main

  implicit none

  integer, parameter :: ni = 96
  integer :: i, ialt
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
  end do
  close(11)

  open(21, file = "./density/COi.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), 1.0d+3
  end do
  close(21)
  open(21, file = "./density/CO.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), 1.0d+3
  end do
  close(21)
  open(21, file = "./density/Ci.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), 1.0d+3
  end do
  close(21)
  open(21, file = "./density/Oi2D.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), 1.0d+3
  end do
  close(21)
  open(21, file = "./density/Oi2P.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), 1.0d+3
  end do
  close(21)
  open(21, file = "./density/N2i.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), 1.0d+3
  end do
  close(21)
  open(21, file = "./density/Ni.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), 1.0d+3
  end do
  close(21)
  open(21, file = "./density/N.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), 1.0d+3
  end do
  close(21)
  open(21, file = "./density/NOi.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), 1.0d+3
  end do
  close(21)
  open(21, file = "./density/HCOi.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), 1.0d+3
  end do
  close(21)
  open(21, file = "./density/N2i.dat", form ="formatted", status = "unknown")
  do i = 1, ni
    write(21, *) altitude(i), 1.0d+3
  end do
  close(21)

end program main
