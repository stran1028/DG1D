module overset
  use code_types
  use pde
contains
  !
  subroutine connect(msh1,msh2)
    implicit none
    type(mesh), intent(inout) :: msh1
    type(mesh), intent(in)    :: msh2
    !
    integer :: i,j
    real*8  :: dx
    !
    ! brute force now 
    !
    do i=1,msh1%nnodes
       do j=1,msh2%nelem
          if (msh2%xe(2,j) .ge. msh1%x(i) .and. &
               msh2%xe(1,j) .le. msh1%x(i)) then
             dx=msh2%xe(2,j)-msh2%xe(1,j)
             if (msh1%nres(i) > dx) then ! blank out cells where m2 is finer
                msh1%iblank(i)=-1
             endif
          endif
       enddo
    enddo
  end subroutine connect
  !
  ! msh2 is considered to be the finer mesh here
  !
  subroutine fixOverlap(msh1,msh2)
    implicit none
    type(mesh), intent(inout) :: msh1
    type(mesh), intent(in) :: msh2
    !
    integer :: i,j,ib1,ib2,ip
    real*8 :: x1,y1,y2
    real*8 :: TOL=1e-8
    integer, parameter :: npass=2
    !
    do ip=1,npass
    m1loop: do i=1,msh1%nnodes
       if (msh1%iblank(i) .ne. 1) cycle m1loop  ! skip blanked nodes
       x1=msh1%x(i)                             ! grab m1 node x coord
       m2loop:do j=1,msh2%nelem
          ib1=msh2%iblank(msh2%e2n(1,j))        ! grab left m2 node
          ib2=msh2%iblank(msh2%e2n(2,j))        ! grab right m2 node
          if (ib1*ib2 > 0 ) cycle m2loop ! skip if either nodes are iblanked 
          y1=msh2%xe(1,j)       ! grab left m2 x coord
          y2=msh2%xe(2,j)       ! grab right m2 x coord
          if ((x1-y1)*(x1-y2) < TOL) then       ! blank m1 node if they're close to overlapping m2 nodes
             msh1%iblank(i)=-1
             cycle m1loop               ! skip the rest of the m2 nodes
          endif
       end do m2loop
    enddo m1loop
    enddo
  end subroutine fixOverlap       
  !
  subroutine findIncompleteElements(msh,elemInfo,nincomp)
    !
    implicit none
    !
    type(mesh), intent(in) :: msh
    integer, intent(out) :: nincomp
    real*8, allocatable, intent(out) :: elemInfo(:)
    integer :: n1,n2,ib1,ib2,i,j,k,nrows
    !
    nincomp=0
    !
    ! find number of elems with only one blanked node
    do i=1,msh%nelem
       n1=msh%e2n(1,i) ! left elem node number
       n2=msh%e2n(2,i) ! right elem node number
       ib1=msh%iblank(n1)
       ib2=msh%iblank(n2)
       if ( ib1*ib2 .le. 0) then
          nincomp=nincomp+1
       endif
    enddo
    !
    ! store info on the incomplete elems
    ! modified for DG, need both q endpts 
    ! and rhs vector
!    nrows = 3 + 2*msh%nshp
!    allocate(elemInfo(nrows*nincomp))
    allocate(elemInfo(nincomp))
!    k=0
    k = 1
    do i=1,msh%nelem
       n1=msh%e2n(1,i)
       n2=msh%e2n(2,i)
       ib1=msh%iblank(n1)
       ib2=msh%iblank(n2)
       if ( ib1*ib2 .le. 0) then
          elemInfo(k) = i
          k=k+1
       endif
    enddo
    !
  end subroutine findIncompleteElements
  !
  subroutine fixFluxIncompleteElements(mshB,mshA,elemInfo,nincomp)
    use bases

    ! Subtract half of overlap section from mesh A (stored in elemInfo)
    implicit none
    type(mesh), intent(inout) :: mshA,mshB
    integer, intent(in) :: nincomp
    real*8, intent(inout) :: elemInfo(nincomp)
!    real*8, intent(inout) :: elemInfo((3+2*msh%nshp)*nincomp)
    !
    integer :: i,j,k,e,nrows,aa,bb,cc,eid
    real*8 :: x1,x2,f1,f2,y1,y2,qA(mshA%nshp),qB(mshB%nshp)
    real*8 :: xcut(2),xc,lcut,xg,vol,flx,ql,qr
    real*8 :: wtmp(mshA%nshp),dwtmp(mshA%nshp)
    real*8 :: qtmp(mshA%nshp),dqtmp(mshA%nshp),dq

!!! XXX add in volume and flux integral for the overset element here

    !
    ! elemInfo = incomplete elements on mesh A
    ! msh = mesh info of mesh B
    !
     iloop: do i=1,nincomp       ! Loop through incomplete elem of mesh A
       eid = elemInfo(i)
       x1=mshA%xe(1,eid) !elemInfo(nrows*(i-1)+2) 
       x2=mshA%xe(2,eid) !elemInfo(nrows*(i-1)+3)
       qA=mshA%q(1,:,eid) !elemInfo(nrows*(i-1)+4:nrows*(i-1)+3+msh%nshp)! 
       eloop: do j=1,mshB%nelem  ! Loop through all elem of mesh B
          y1=mshB%xe(1,j)
          y2=mshB%xe(2,j)
          if (x1 > y2 .or. x2 < y1) then ! skip if not overlapping
             cycle eloop
          else
	     if (mshB%iblank(mshB%e2n(1,j)) .ne.1 .and. &
                 mshB%iblank(mshB%e2n(2,j)) .ne.1) cycle eloop ! skip if incomplete mesh B elem
             qB=mshB%q(1,:,j)

             if ((x1-y1)*(x1-y2) .le. 0.0) then ! L node of mesh A is inside of mesh B elem
               ! Overlap is between x1 and y2
               ! mshA will remove first half of overlap (from x1 to 0.5*(x1+y2))
               xcut = [x1,0.5d0*(x1+y2)]
               
               ! add intermesh flux from mesh B interior to mesh A L node
               call shapefunction(mshA%nshp,-0.5d0,[-0.5d0,0.5d0],[1d0,1d0],qtmp,dqtmp)
               wtmp = qtmp
               call shapefunction(mshB%nshp,x1,[y1,y2],qB,qtmp,dqtmp)
               qL = sum(qtmp)
               call shapefunction(mshA%nshp,x1,[x1,x2],qA,qtmp,dqtmp)
               qR = sum(qtmp)
               call flux(ql,qr,flx)
               do k = 1,mshA%nshp
                 mshA%rhs(:,k,eid) = mshA%rhs(:,k,eid) + wtmp(k)*flx
               enddo

             elseif ((x2-y1)*(x2-y2) .le. 0.0) then ! R ndoe of mesh A is inside of mesh B elem          
               ! Overlap is between y1 and x2
               ! msh A will remove second half of overlap (from 0.5(y1+x2) to x2
               xcut = [0.5d0*(y1+x2),x2]

               ! add intermesh flux from mesh B interior to mesh A R node
               call shapefunction(mshA%nshp,0.5d0,[-0.5d0,0.5d0],[1d0,1d0],qtmp,dqtmp)
               wtmp = qtmp
               call shapefunction(mshB%nshp,x2,[y1,y2],qB,qtmp,dqtmp)
               qR = sum(qtmp)
               call shapefunction(mshA%nshp,x2,[x1,x2],qA,qtmp,dqtmp)
               qL = sum(qtmp)
               call flux(ql,qr,flx)
               do k = 1,mshA%nshp
                 mshA%rhs(:,k,eid) = mshA%rhs(:,k,eid) - wtmp(k)*flx
               enddo

             endif
             lcut = xcut(2)-xcut(1)
             xc = 0.5*(xcut(1)+xcut(2))   ! center of section to be removed
               
             ! Adjust mass matrix and volume integral
             do aa = 1,mshA%ngauss
               ! get shapefunction from msh A at quad pts of cut section
               xg = mshA%xgauss(aa)*lcut+xc
               call shapefunction(mshA%nshp,xg,[x1,x2],[1d0,1d0],wtmp,dwtmp)
               call shapefunction(mshA%nshp,xg,[x1,x2],qA,qtmp,dqtmp)
               do bb = 1,mshA%nshp
                 ! Volume Integral
                 ! Remove int(w u,x dx) 
                 ! Note that we are removing the original volume integral
                 ! (before integration by parts) since it's easier and 
                 ! mathematically the same 
                 call volint(sum(dqtmp),vol)
                 mshA%rhs(:,bb,eid) = mshA%rhs(:,bb,eid) + wtmp(bb)*vol*mshA%wgauss(aa)*lcut/mshA%dx(eID) ! scale gauss weights by cut length in parent element

               enddo ! nshp
             enddo ! ngauss
             cycle iloop
          endif
       enddo eloop
    enddo iloop
  end subroutine fixFluxIncompleteElements
  !
  subroutine fixMassIncompleteElements(mshB,mshA,elemInfo,nincomp)
    use bases

    ! Subtract half of overlap section from mesh A (stored in elemInfo)
    implicit none
    type(mesh), intent(inout) :: mshA,mshB
    integer, intent(in) :: nincomp
    real*8, intent(inout) :: elemInfo(nincomp)
!    real*8, intent(inout) :: elemInfo((3+2*msh%nshp)*nincomp)
    !
    integer :: i,j,k,e,nrows,aa,bb,cc,eid,index1
    real*8 :: x1,x2,f1,f2,y1,y2,qA(mshA%nshp),qB(mshB%nshp)
    real*8 :: xcut(2),xc,lcut,xg
    real*8 :: wtmp(mshA%nshp),dwtmp(mshA%nshp)
    !
    ! elemInfo = incomplete elements on mesh A
    ! msh = mesh info of mesh B
    !
     iloop: do i=1,nincomp       ! Loop through incomplete elem of mesh A
       eid = elemInfo(i)
       x1=mshA%xe(1,eid) !elemInfo(nrows*(i-1)+2) 
       x2=mshA%xe(2,eid) !elemInfo(nrows*(i-1)+3)
       qA=mshA%q(1,:,eid) !elemInfo(nrows*(i-1)+4:nrows*(i-1)+3+msh%nshp)! 
       eloop: do j=1,mshB%nelem  ! Loop through all elem of mesh B
          y1=mshB%xe(1,j)
          y2=mshB%xe(2,j)
          if (x1 > y2 .or. x2 < y1) then ! skip if not overlapping
             cycle eloop
          else
	     if (mshB%iblank(mshB%e2n(1,j)) .ne.1 .and. &
                 mshB%iblank(mshB%e2n(2,j)) .ne.1) cycle eloop ! skip if incomplete mesh B elem
             qB=mshB%q(1,:,j)

             if ((x1-y1)*(x1-y2) .le. 0.0) then ! L node of mesh A is inside of mesh B elem
               ! Overlap is between x1 and y2
               ! mshA will remove first half of overlap (from x1 to 0.5*(x1+y2))
               xcut = [x1,0.5d0*(x1+y2)]

             elseif ((x2-y1)*(x2-y2) .le. 0.0) then ! R ndoe of mesh A is inside of mesh B elem          
               ! Overlap is between y1 and x2
               ! msh A will remove second half of overlap (from 0.5(y1+x2) to x2
               xcut = [0.5d0*(y1+x2),x2]

             endif
             lcut = xcut(2)-xcut(1)
             xc = 0.5*(xcut(1)+xcut(2))   ! center of section to be removed
               
             ! Adjust mass matrix and volume integral
             do aa = 1,mshA%ngauss
               ! get shapefunction from msh A at quad pts of cut section
               xg = mshA%xgauss(aa)*lcut+xc
               call shapefunction(mshA%nshp,xg,[x1,x2],[1d0,1d0],wtmp,dwtmp)
               do bb = 1,mshA%nshp
                 do cc = 1,mshA%nshp
                    index1 = (bb-1)*mshA%nshp+cc
                    ! Fix mass matrix
                    mshA%mass(:,index1,eid) = mshA%mass(:,index1,eid) - wtmp(bb)*wtmp(cc)*mshA%wgauss(aa)*lcut
                 enddo ! nshp

               enddo ! nshp
             enddo ! ngauss

             cycle iloop
          endif
       enddo eloop
    enddo iloop
  end subroutine fixMassIncompleteElements
  !
end module overset
          
  
