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
    ! Verify this still makes sense with DG indices XXX
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
!          elemInfo(nrows*k+1) = i
!          elemInfo(nrows*k+2) = msh%xe(1,i)
!          elemInfo(nrows*k+3) = msh%xe(2,i)
!          do j=4,3+msh%nshp
!            elemInfo(nrows*k+j) = msh%q(1,j-3,i)
!          enddo
!          do j=4+msh%nshp,3+2*msh%nshp
!            elemInfo(nrows*k+j) = msh%rhs(1,j-3-msh%nshp,i)
!          enddo
          k=k+1
       endif
    enddo
    !
  end subroutine findIncompleteElements
  !
  subroutine fixFluxIncompleteElements(mshB,mshA,elemInfo,nincomp)
    ! Subtract half of overlap section from mesh A (stored in elemInfo)
    implicit none
    type(mesh), intent(inout) :: mshA,mshB
    integer, intent(in) :: nincomp
    real*8, intent(inout) :: elemInfo(nincomp)
!    real*8, intent(inout) :: elemInfo((3+2*msh%nshp)*nincomp)
    !
    integer :: i,j,e,nrows,aa,bb,cc,eid
    real*8 :: x1,x2,f1,f2,y1,y2,qA(mshA%nshp),qB(mshB%nshp)
    real*8 :: dx,xs,qtmp(mshA%nshp),dqtmp(mshA%nshp),qs,vol,qval
    real*8 :: temp(2),xg

!! Do I need all of these elemInfo arrays? Don't I just need i_incomp? It's
!! confusing this way
    !
    ! elemInfo = incomplete elements on mesh A
    ! msh = mesh info of mesh B
    !
!    nrows = 3 + 2*msh%nshp              
     iloop: do i=1,nincomp       ! Loop through incomplete elem of mesh A
       write(*,*) 'ST DEBUG eINFO = ',elemInfo
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
               ! Get the middle point of the overlap
               ! Overlap is between x1 and y2
               xs = 0.5*(x1+y2) ! overlap split point
               dx = abs(x1-xs)  ! subtracting half of the oberlap
               call shapefunction(mshA%nshp,xs,[x1, x2],qA,qtmp,dqtmp)
               qs = SUM(qtmp)              

               ! Adjust mass matrix
               do aa = 1,mshA%ngauss
                 xg = (mshA%xgauss(aa)+0.5)*dx + xs
                 call shapefunction(mshA%nshp,xg,[xs,x2],[1,1],qtmp,dqtmp)
                 do bb=1,mshA%nshp
                 do cc=1,mshA%nshp
                   mshA%mass(1,bb,cc,eid) = mshA%mass(1,bb,cc,eid) - qtmp(bb)*qtmp(cc)*dx*mshA%wgauss(aa)
                 enddo
                 enddo
               enddo

               ! Calculate the volume integral in overlap element
               temp = 0d0
               do aa = 1,mshA%ngauss
                 qtmp = 0.0
                 ! wrong! xgauss and x1,xs have to match up. need to transform
                 ! xgauss or vice versa
                 xg = (mshA%xgauss(aa)+0.5)*dx + x1
                 call shapefunction(mshA%nshp,xg,[x1,xs],[qA(1),qs],qtmp,dqtmp)
                 qval = SUM(qtmp)
                 do bb = 1,mshA%nshp
                    call volint(qval,vol)
                    vol = vol*mshA%dshp(aa,bb)*dx*mshA%wgauss(aa)
                    temp(bb) = temp(bb) + vol
                    mshA%rhs(1,:,eid) = mshA%rhs(1,:,eid) + vol 
                    !elemInfo(3+mshA%nshp+bb) = elemInfo(3+mshA%nshp+bb) + vol 
                 enddo
               enddo

               ! Calculate flux vector in overlap element
               ! Subtract from left side of mesh A element
               call flux(qA(1),qA(1),f1)
               call flux(qs,qs,f2)
               mshA%rhs(1,1,eid) = mshA%rhs(1,1,eid) - (f2-f1)
               !elemInfo(4+mshA%nshp) = elemInfo(4+msh%nshp) - (f2-f1) 
                 write(*,*) ' L Side'
                 write(*,*) '  i,j =',i,j
                 write(*,*) '  xs,dx = ',xs,dx
                 write(*,*) '  x1,x2 = ',x1,x2
                 write(*,*) '  y1,y2 = ',y1,y2
                 write(*,*) '  qval = ',qval    ! wrong
                 write(*,*) '  qA = ',qA(1),qA(2)
                 write(*,*) '  qB = ',qB(1),qB(2)
                 write(*,*) '  vol = ',temp(1),temp(2)
                 write(*,*) '  f1,f2 = ', f1,f2
                 write(*,*) ' '

               cycle iloop
             endif
             if ((x2-y1)*(x2-y2) .le. 0.0) then ! R node of mesh A is inside of mesh B elem
               ! overlap is between y1 and x2
               xs = 0.5*(y1+x2) ! overlap split point
               dx = abs(y1-xs)  ! subtracting half of the oberlap
               call shapefunction(mshA%nshp,xs,[x1, x2],qA,qtmp,dqtmp)
               qs = SUM(qtmp)             

               ! Adjust mass matrix
               do aa = 1,mshA%ngauss
                 xg = (mshA%xgauss(aa)+0.5)*dx + xs
                 call shapefunction(mshA%nshp,xg,[xs,x2],[1,1],qtmp,dqtmp)
                 do bb=1,mshA%nshp
                 do cc=1,mshA%nshp
                   mshA%mass(1,bb,cc,eid) = mshA%mass(1,bb,cc,eid) - qtmp(bb)*qtmp(cc)*dx*mshA%wgauss(aa)
                 enddo
                 enddo
               enddo

               ! Calculate the volume integral in overlap elemen
               temp = 0d0
               do aa = 1,mshA%ngauss
                 xg = (mshA%xgauss(aa)+0.5)*dx + xs
                 call shapefunction(mshA%nshp,xg,[xs,x2],[qs,qA(mshA%nshp)],qtmp,dqtmp)
                 qval = SUM(qtmp)
                 do bb = 1,mshA%nshp
                    call volint(qval,vol)
                    vol = vol*mshA%dshp(aa,bb)*dx*mshA%wgauss(aa)
                    temp(bb) = temp(bb)+vol
                    mshA%rhs(1,:,eid) = mshA%rhs(1,:,eid) + vol
!                    elemInfo(3+mshA%nshp+bb) = elemInfo(3+mshA%nshp+bb) + vol
                 enddo
               enddo

               ! Calculate flux vector in overlap element
               ! Subtract from right side of mesh A element
               call flux(qA(mshA%nshp),qA(mshA%nshp),f1)
               call flux(qs,qs,f2)
               mshA%rhs(1,mshA%nshp,eid) = mshA%rhs(1,mshA%nshp,eid) - (f2-f1)
!               elemInfo(3+2*mshA%nshp) = elemInfo(3+2*mshA%nshp) - (f2-f1)
                 write(*,*) ' R Side'
                 write(*,*) '  i,j =',i,j
                 write(*,*) '  xs,dx = ',xs,dx
                 write(*,*) '  x1,x2 = ',x1,x2
                 write(*,*) '  y1,y2 = ',y1,y2
                 write(*,*) '  qtmp = ',qtmp
                 write(*,*) '  qA = ',qA(1),qA(2)
                 write(*,*) '  qB = ',qB(1),qB(2)
                 write(*,*) '  vol = ',temp(1),temp(2)
                 write(*,*) '  f1,f2 = ', f1,f2
                 write(*,*) ' '

               cycle iloop
             endif
          endif
       enddo eloop
    enddo iloop
  end subroutine fixFluxIncompleteElements
  !
  subroutine setRHS(msh,elemInfo,nincomp)
    implicit none
    type(mesh), intent(inout) :: msh
    integer :: nincomp
    real*8, intent(in) :: elemInfo((3+2*msh%nshp)*nincomp)
    integer :: i,e,nrows
    real*8 :: x1,x2,res(msh%nshp)
    integer, save :: iout=0
    !
    nrows = 3 + 2*msh%nshp
    ! Replace the fluxes for the incomplete elements
    ! Need to fix for DG XXX
    do i=1,nincomp
       e=nint(elemInfo(nrows*(i-1)+1))
       x1=elemInfo(nrows*(i-1)+2)
       x2=elemInfo(nrows*(i-1)+3)
       !write(6,*) 'e,x1,x2=',e,x1,x2,msh%rhs(1,1,e)
       !write(6,*) 'msh%xe=',msh%xe(:,e)
       !write(6,*) 'r1:',msh%rhs(1,1,e)
       msh%xe(1,e)=x1
       msh%xe(2,e)=x2
       msh%rhs(1,:,e)=elemInfo(4+msh%nshp:3+2*msh%nshp)
       !write(6,*) 'r2:',msh%rhs(1,1,e)
       !write(6,*) 'e,x1,x2=',e,x1,x2,msh%rhs(1,1,e)
    enddo
    !
    !iout=iout+1
    !do i=1,msh%nelem
    !   write(40+iout,*) msh%rhs(1,1,i),msh%iblank(msh%e2n(1,i)),&
    !        msh%iblank(msh%e2n(2,i))
    !enddo
    !
  end subroutine setRHS
  !
end module overset
          
  
