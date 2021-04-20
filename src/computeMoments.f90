subroutine computeMoments(msh,moments)
 !
  use bases
  use code_types
 !
 implicit none
 type(mesh), intent(inout) :: msh
 real*8, intent(inout) :: moments(2)
 integer :: i,j,k
 real*8 :: dx,xc,qvals(msh%nshp),dqvals(msh%nshp)
 real*8 ::xlen
 !
 moments=0d0
 xlen=0d0
 do i=1,msh%nelem
   if (msh%iblank(1,i) .ne. 1 .and. &
       msh%iblank(2,i) .ne. 1) cycle
   xc=(msh%xe(2,i)+msh%xe(1,i))*0.5

    do j = 1,msh%ngauss
      qvals = 0d0
      call shapefunction(msh%nshp,msh%xgauss(j),[-0.5d0,0.5d0],msh%q(1,:,i),qvals,dqvals)
      moments(1) = moments(1) + msh%wgauss(j)*sum(qvals)*msh%dx(i) 
      moments(2) = moments(2) + msh%wgauss(j)*sum(qvals)*msh%dx(i)*xc
    enddo

 enddo
! moments(2)=xlen
write(*,*) ' '
write(*,*) 'Moments: = ',moments
write(*,*) ' '
end subroutine computeMoments
