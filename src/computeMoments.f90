subroutine computeMoments(msh,moments,error,nincomp,elemInfo)
 !
  use pde
  use bases
  use code_types
 !
 implicit none
 integer,intent(in) :: nincomp
 real*8,intent(in) :: elemInfo(3,nincomp)
 type(mesh), intent(inout) :: msh
 real*8, intent(inout) :: moments(2)
 integer :: i,j,k,ndof,eid
 real*8 :: dx,xc,qvals(msh%nshp),dqvals(msh%nshp)
 real*8 ::xlen,xg,exact(1),error(1),tmp
 !
 moments=0d0
 xlen=0d0
 ndof = 0
 error = 0d0
 ! Interior points
 do i=1,msh%nelem
! if (msh%iblank(1,i) .ne. 1 .and. &
!       msh%iblank(2,i) .ne. 1) cycle
   if(minval(msh%iblank(:,i)).lt.0) cycle

     xc=(msh%xe(2,i)+msh%xe(1,i))*0.5
     ndof = ndof + msh%nshp
     tmp = 0d0
     do j = 1,msh%ngauss
       qvals = 0d0
       call shapefunction(msh%nshp,msh%xgauss(j),[-0.5d0,0.5d0],msh%q(1,:,i),qvals,dqvals)
       tmp = tmp+msh%wgauss(j)*sum(qvals)*msh%dx(i)
       moments(1) = moments(1) + msh%wgauss(j)*sum(qvals)*msh%dx(i) 
       moments(2) = moments(2) + msh%wgauss(j)*sum(qvals)*msh%dx(i)*xc
 
       xg = msh%xgauss(j)*dx + xc
       call initq(xg,exact)
       error = error + msh%wgauss(j)*(exact-sum(qvals))**2d0*msh%dx(i)
     enddo
 enddo
 ! Fringe elements
 do i = 1,nincomp
     ! Get the modified element size
     eid = elemInfo(1,i)
     dx = elemInfo(3,i)-elemInfo(2,i)
     xc = 0.5d0*(elemInfo(3,i)+elemInfo(2,i))

     ! adjust gauss points to modified element
     ndof = ndof + msh%nshp
     tmp = 0d0
     do j = 1,msh%ngauss
       xg = msh%xgauss(j)*dx + xc
       qvals = 0d0
       call shapefunction(msh%nshp,xg,[msh%xe(1,eid),msh%xe(2,eid)],msh%q(1,:,eid),qvals,dqvals)
       tmp = tmp+msh%wgauss(j)*sum(qvals)*msh%dx(i)
       moments(1) = moments(1) + msh%wgauss(j)*sum(qvals)*dx
       moments(2) = moments(2) + msh%wgauss(j)*sum(qvals)*dx*xc

       call initq(xg,exact)
       error = error + msh%wgauss(j)*(exact-sum(qvals))**2d0*msh%dx(i)
     enddo
 enddo
 
! moments(2)=xlen
!write(*,*) ' '
!write(*,*) 'ndof, Moments: = ',ndof,moments
!write(*,*) 'L2 error: ',error
!write(*,*) ' '
end subroutine computeMoments
