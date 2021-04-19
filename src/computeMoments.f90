subroutine computeMoments(msh,moments)
 !
 use code_types
 !
 implicit none
 type(mesh), intent(inout) :: msh
 real*8, intent(inout) :: moments(2)
 integer :: i
 real*8 :: dx,xc
 real*8 ::xlen
 !
 moments=0d0
 xlen=0d0
 do i=1,msh%nelem
   if (msh%iblank(1,i) .ne. 1 .and. &
       msh%iblank(2,i) .ne. 1) cycle
   dx=msh%xe(2,i)-msh%xe(1,i)
   xc=(msh%xe(2,i)+msh%xe(1,i))*0.5
   moments(1)=moments(1)+msh%q(1,1,i)*dx
   moments(2)=moments(2)+msh%q(1,1,i)*dx*xc
   xlen=xlen+dx
 enddo
 moments(2)=xlen
 !write(6,*) xlen,moments(1)
end subroutine computeMoments
