subroutine shapefunction(nshp,x,xlim,q,qout,dqout)
    implicit none
  !
  integer, intent(in) :: nshp
  real*8, intent(in) :: xlim(2),x,q(nshp)
  real*8, intent(inout) :: qout(nshp),dqout(nshp)
  !
  integer:: i,j,k
  real*8 :: dx,xloc
  !
  ! Transform element to go from [-.5 .5]
  dx = xlim(2)-xlim(1)
  xloc = (x-xlim(1))/dx - 0.5
  !
  if(nshp.eq.2) then ! 1st order Lagrange elements
    qout(1)  = q(1)*(0.5-xloc)
    dqout(1) = -1d0*q(1)
    qout(2)  = q(2)*(0.5+xloc)
    dqout(2) = 1d0*q(2)
  endif 

end subroutine shapefunction
