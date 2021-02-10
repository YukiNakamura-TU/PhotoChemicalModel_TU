module p__airglow
  use v__tdec,        only : spl_, var_, grd_, cst_, xct_, flx_

  implicit none
  integer(4), parameter :: sp = 4, dp = 8

  private
  public :: p__airglow__exe

contains


  subroutine p__airglow__exe(spl, cst, grd, & ! in
      &                      var            ) ! inout
    type(cst_),           intent(in)    :: cst
    type(spl_),           intent(in)    :: spl
    type(grd_),           intent(in)    :: grd
    type(var_),           intent(inout) :: var

  end subroutine p__airglow__exe


end module p__airglow
