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
 integer :: i,j,k,ndof,eid,pid
 real*8 :: dxmod,xc,qvals(msh%nshp),dqvals(msh%nshp)
 real*8 :: q0vals(msh%nshp),dq0vals(msh%nshp)
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
       moments(1) = moments(1) + msh%wgauss(j)*sum(qvals)*msh%dx(i) 
       moments(2) = moments(2) + msh%wgauss(j)*sum(qvals)*msh%dx(i)*xc
 
       ! exact solution held in q0
       q0vals = 0d0
       call shapefunction(msh%nshp,msh%xgauss(j),[-0.5d0,0.5d0],msh%q0(1,:,i),q0vals,dq0vals)
       error = error + msh%wgauss(j)*(sum(q0vals)-sum(qvals))**2d0*msh%dx(i)
     enddo
 enddo
 ! Fringe elements
 do i = 1,nincomp
     ! Get the modified element size
     eid = elemInfo(1,i)
     pid = msh%parent(eid)
     dxmod = elemInfo(3,i)-elemInfo(2,i)
     xc = 0.5d0*(elemInfo(3,i)+elemInfo(2,i))

     ! adjust gauss points to modified element
     ndof = ndof + msh%nshp
     tmp = 0d0
     do j = 1,msh%ngauss
       xg = msh%xgauss(j)*dxmod + xc
       qvals = 0d0
       call shapefunction(msh%nshp,xg,[msh%xe(1,pid),msh%xe(2,pid)],msh%q(1,:,pid),qvals,dqvals)
       tmp = tmp+msh%wgauss(j)*sum(qvals)*msh%dx(i)
       moments(1) = moments(1) + msh%wgauss(j)*sum(qvals)*dxmod
       moments(2) = moments(2) + msh%wgauss(j)*sum(qvals)*dxmod*xc

       ! exact solution held in q0
       q0vals = 0d0
       call shapefunction(msh%nshp,xg,[msh%xe(1,pid),msh%xe(2,pid)],msh%q0(1,:,pid),q0vals,dq0vals)
       error = error + msh%wgauss(j)*(sum(q0vals)-sum(qvals))**2d0*dxmod
     enddo
 enddo
 
end subroutine computeMoments
