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
    integer :: i,j,k,index1
    real*8  :: dx,res
    !
    ! brute force now 
    !
    do i = 1,msh1%nelem
    do j = 1,msh2%nelem
      dx = msh2%xe(2,j)-msh2%xe(1,j)

      !If both vtx of mesh 1 in same element of mesh 2
      if((msh2%xe(2,j).ge.(msh1%xe(2,i))).and. &
         (msh2%xe(1,j).le.(msh1%xe(2,i))).and. &
         (msh2%xe(2,j).ge.(msh1%xe(1,i))).and. &
         (msh2%xe(1,j).le.(msh1%xe(1,i)))) then 
         index1 = 2*(i-1)

         ! If 1 node needs to be blanked, both do
         res = maxval([msh1%nres(index1+1),msh1%nres(index1+2)])
         if(res.gt.dx) msh1%iblank(:,i) = -1
      else ! only 1 vtx overlapping in element of mesh 2
         do k = 1,2
           if (msh2%xe(2,j) .ge. msh1%xe(k,i) .and. &
                msh2%xe(1,j) .le. msh1%xe(k,i)) then


              index1 = 2*(i-1)+k
              if(msh1%nres(index1).gt.dx) msh1%iblank(k,i)=-1
           endif
         enddo
      endif
    enddo ! mesh 2 elem
    enddo ! mesh 1 elem

    ! Go back and fix any mistakes
    do i = 2,msh1%nelem-1
      ! If previous element was all blanked, 
      ! L node of this elem needs to be blanked
      ! leaving only a partial overlap
      if(msh1%iblank(2,i-1).eq.-1 ) then   
        msh1%iblank(1,i)=-1
      endif

      ! If next element was all blanked, 
      ! R node of this elem needs to be blanked
      ! leaving only a partial overlap
      if(msh1%iblank(1,i+1).eq.-1 ) then   
        msh1%iblank(2,i)=-1
      endif
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
    integer :: i,j,ib1,ib2,ip,aa,bb
    real*8 :: x1,y1,y2
    real*8 :: TOL=1e-8
    integer, parameter :: npass=2
    !
    do ip=1,npass
!    m1loop: do i=1,msh1%nnodes
    m1loop: do aa=1,msh1%nelem
       do bb = 1,2
         if (msh1%iblank(bb,aa) .ne. 1) cycle m1loop  ! skip blanked nodes
         x1=msh1%xe(bb,aa)                             ! grab m1 node x coord
         m2loop:do j=1,msh2%nelem
            ib1=msh2%iblank(1,j)        ! grab left m2 node
            ib2=msh2%iblank(2,j)        ! grab right m2 node
            if (ib1*ib2 > 0 ) cycle m2loop ! skip if either nodes are iblanked 
            y1=msh2%xe(1,j)       ! grab left m2 x coord
            y2=msh2%xe(2,j)       ! grab right m2 x coord
            if ((x1-y1)*(x1-y2) < TOL) then       ! blank m1 node if they're close to overlapping m2 nodes
               msh1%iblank(bb,aa)=-1
               cycle m1loop               ! skip the rest of the m2 nodes
            endif
         end do m2loop
      enddo
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
    real*8, allocatable, intent(out) :: elemInfo(:,:)
    integer :: n1,n2,ib1,ib2,i,j,k,nrows
    !
    nincomp=0
    !
    ! find number of elems with only one blanked node
    do i=1,msh%nelem
       ib1=msh%iblank(1,i)
       ib2=msh%iblank(2,i)
       if ( ib1*ib2 .le. 0) then
          nincomp=nincomp+1
       endif
    enddo
    !
    ! store info on the incomplete elems
    ! modified for DG, need both q endpts 
    ! and rhs vector
    allocate(elemInfo(3,nincomp))
    k = 1
    do i=1,msh%nelem
       ib1=msh%iblank(1,i)
       ib2=msh%iblank(2,i)
       if ( ib1*ib2 .le. 0) then
          elemInfo(1,k) = i
          k=k+1
       endif
    enddo
    !
  end subroutine findIncompleteElements
  !
  subroutine fixFluxIncompleteElements(mshB,mshA,elemInfo,nincomp,consoverset,foverlap,isupg,dt)
    use bases

    ! Subtract half of overlap section from mesh A (stored in elemInfo)
    implicit none
    type(mesh), intent(inout) :: mshA,mshB
    integer, intent(in) :: nincomp,consoverset,isupg
    real*8, intent(inout) :: elemInfo(3,nincomp),foverlap
    real*8, intent(in) :: dt
    !
    integer :: i,j,k,e,nrows,aa,bb,cc,eid,neigh
    real*8 :: x1,x2,f1,f2,y1,y2,qA(mshA%nshp),qB(mshB%nshp)
    real*8 :: xrem(2),xcut(2),xc,lcut,xg,vol,flx,qL,qR,fact,xfac
    real*8 :: wtmp(mshA%nshp),dwtmp(mshA%nshp),tmp
    real*8 ::qtmp(mshA%nshp),dqtmp(mshA%nshp),dq,dvol(mshA%nshp),dflx(mshA%nshp)
    real*8 :: qval,dqval,dudt,resid,tau
    !
    ! elemInfo = incomplete elements on mesh A
    ! msh = mesh info of mesh B
    !
     iloop: do i=1,nincomp       ! Loop through incomplete elem of mesh A
       eid = elemInfo(1,i)
       x1=mshA%xe(1,eid) 
       x2=mshA%xe(2,eid) 
       qA=mshA%q(1,:,eid) 
       eloop: do j=1,mshB%nelem  ! Loop through all elem of mesh B
          y1=mshB%xe(1,j)
          y2=mshB%xe(2,j)
          if (x1 > y2 .or. x2 < y1) then ! skip if not overlapping
             cycle eloop
          else
	     if (mshB%iblank(1,j) .ne.1 .and. &
                 mshB%iblank(2,j) .ne.1) cycle eloop ! skip if element blanked
             qB=mshB%q(1,:,j)

             if(mshA%dx(eid).lt.mshB%dx(j)) then ! A is fine mesh
               xfac = foverlap
             else ! A is coarse mesh
               xfac = 1d0-foverlap
             endif

             if ((x1-y1)*(x1-y2) .le. 0.0) then ! L node of mesh A is inside of mesh B elem
               ! Full overlap is between x1 and y2
               ! mshA will remove first section of overlap 
               if(consoverset.eq.1) then 
                 xcut = [x1,x1+xfac*(y2-x1)]
               else
                 xcut = [x1,x1]
               endif
               xrem = elemInfo(2:3,i)

               ! add intermesh flux from mesh B interior to mesh A L node
               wtmp = 1d0
               call shapefunction(mshA%nshp,xcut(2),[x1,x2],wtmp,wtmp,dwtmp)
               call shapefunction(mshB%nshp,xcut(2),[y1,y2],qB,qtmp,dqtmp)
               qL = sum(qtmp)
               call shapefunction(mshA%nshp,xcut(2),[x1,x2],qA,qtmp,dqtmp)
               qR = sum(qtmp)
               call flux(qL,qR,flx)
               do k = 1,mshA%nshp
                 mshA%rhs(:,k,eid) = mshA%rhs(:,k,eid) + wtmp(k)*flx
               enddo

               ! Handle mesh A R node flux
               wtmp = 1d0
               call shapefunction(mshA%nshp,x2,[x1,x2],wtmp,wtmp,dwtmp)
               neigh = mshA%face(2,eid)
               call shapefunction(mshA%nshp,x2,[x1,x2],qA,qtmp,dqtmp)
               qL = sum(qtmp)
               call shapefunction(mshA%nshp,x2,mshA%xe(:,neigh),mshA%q(1,:,neigh),qtmp,dqtmp)
               qR = sum(qtmp)
               call flux(qL,qR,flx)
               do k = 1,mshA%nshp
                 mshA%rhs(:,k,eid) = mshA%rhs(:,k,eid) - wtmp(k)*flx
               enddo

             elseif ((x2-y1)*(x2-y2) .le. 0.0) then ! R node of mesh A is inside of mesh B elem          
               ! Overlap is between y1 and x2
               ! msh A will remove second half of overlap (from 0.5(y1+x2) to x2
               if(consoverset.eq.1) then 
                 xcut = [x2-xfac*(x2-y1),x2]
               else
                 xcut = [x2,x2]
               endif
               xrem = elemInfo(2:3,i)

               ! Handle mesh A L node flux 
               wtmp = 1d0
               call shapefunction(mshA%nshp,x1,[x1,x2],wtmp,wtmp,dwtmp)
               neigh = mshA%face(1,eid)
               call shapefunction(mshA%nshp,x1,mshA%xe(:,neigh),mshA%q(1,:,neigh),qtmp,dqtmp)
               qL = sum(qtmp)
               call shapefunction(mshA%nshp,x1,[x1,x2],qA,qtmp,dqtmp)
               qR = sum(qtmp)
               call flux(qL,qR,flx)
               do k = 1,mshA%nshp
                 mshA%rhs(:,k,eid) = mshA%rhs(:,k,eid) + wtmp(k)*flx
               enddo

               ! add intermesh flux from mesh B interior to mesh A R node
               wtmp = 1d0
               call shapefunction(mshA%nshp,xcut(1),[x1,x2],wtmp,wtmp,dwtmp)
               call shapefunction(mshB%nshp,xcut(1),[y1,y2],qB,qtmp,dqtmp)
               qR = sum(qtmp)
               call shapefunction(mshA%nshp,xcut(1),[x1,x2],qA,qtmp,dqtmp)
               qL = sum(qtmp)
               call flux(ql,qr,flx)
               do k = 1,mshA%nshp
                 mshA%rhs(:,k,eid) = mshA%rhs(:,k,eid) - wtmp(k)*flx
               enddo

             endif
             fact = (xrem(2)-xrem(1))/(x2-x1)
               
             dvol = 0d0
             ! Compute volume integral over partial element 
             do aa = 1,mshA%ngauss
               ! get shapefunction from msh A at quad pts of remaining element
               xg = mshA%xgauss(aa)*(xrem(2)-xrem(1))+0.5d0*(xrem(2)+xrem(1))
               wtmp = 1d0
               call shapefunction(mshA%nshp,xg,[x1,x2],wtmp,wtmp,dwtmp)
               call shapefunction(mshA%nshp,xg,[x1,x2],qA,qtmp,dqtmp)
               qval = sum(qtmp)
               dqval = sum(dqtmp)/mshA%dx(eid)
               call volint(qval,vol)
               do bb = 1,mshA%nshp
                 ! Volume Integral
                 mshA%rhs(:,bb,eid) = mshA%rhs(:,bb,eid) + dwtmp(bb)*vol*(mshA%wgauss(aa)*fact) ! scale gauss weights by length of remaining element parent

                 dvol(bb) = dvol(bb) + dwtmp(bb)*vol*(mshA%wgauss(aa)*fact)

                 ! SUPG Terms
                 if(isupg.eq.1) then
                   ! get residual
                   call shapefunction(mshA%nshp,xg,[x1,x2],mshA%qold(1,:,eid),qtmp,dqtmp)
                   dudt = (qval-sum(qtmp))/dt
                   if (index(pde_descriptor,'linear_advection') > 0 ) then
                     resid = dudt + a*dqval ! du/dt + a du/dx
                     tau = sqrt(a*a/mshA%dxcut(eid)/mshA%dxcut(eid) + 4d0/dt/dt)
                     tau = a/tau
                   else if (index(pde_descriptor,'burgers') > 0) then
                     resid = dudt + qval*dqval ! du/dt + u du/dx
                     tau = sqrt(qval*qval/mshA%dxcut(eid)/mshA%dxcut(eid) + 4d0/dt/dt)
                     tau = qval/tau
                   endif
                   if(qval.gt.1.5d0) then 
                   write(*,*) ' '
                   write(*,*) ' eid = ',eid
                   write(*,*) '   u,uold,dt = ',qval, sum(qtmp),dt
                   write(*,*) '   dudt, dqval =',dudt, dqval
                   write(*,*) '   resid = ',resid
                   write(*,*) '   tau = ',tau
                   write(*,*) '   supg = ',tau*resid
                   endif
                   mshA%rhs(:,bb,eid) = mshA%rhs(:,bb,eid) - dwtmp(bb)*tau*resid*(mshA%wgauss(aa)*fact)
                 endif ! supg
               enddo ! nshp
             enddo ! ngauss
             cycle iloop
          endif
       enddo eloop
    enddo iloop
  end subroutine fixFluxIncompleteElements
  !
  subroutine fixMassIncompleteElements(mshB,mshA,elemInfo,nincomp,consoverset,foverlap)
    use bases

    ! Subtract half of overlap section from mesh A (stored in elemInfo)
    implicit none
    type(mesh), intent(inout) :: mshA,mshB
    integer, intent(in) :: nincomp,consoverset
    real*8, intent(inout) :: elemInfo(3,nincomp),foverlap
    !
    integer :: i,j,k,e,nrows,aa,bb,cc,eid,index1
    real*8 :: x1,x2,f1,f2,y1,y2,qA(mshA%nshp),qB(mshB%nshp)
    real*8 :: xcut(2),xc,lcut,xg,xfac
    real*8 :: wtmp(mshA%nshp),dwtmp(mshA%nshp)
    !
    ! elemInfo = incomplete elements on mesh A
    ! msh = mesh info of mesh B
    !
     iloop: do i=1,nincomp       ! Loop through incomplete elem of mesh A
       eid = elemInfo(1,i)
       x1=mshA%xe(1,eid) 
       x2=mshA%xe(2,eid) 
       qA=mshA%q(1,:,eid) 
       eloop: do j=1,mshB%nelem  ! Loop through all elem of mesh B
          y1=mshB%xe(1,j)
          y2=mshB%xe(2,j)
          if (x1 > y2 .or. x2 < y1) then ! skip if not overlapping
             cycle eloop
          else
	     if (mshB%iblank(1,j) .ne.1 .and. &
                 mshB%iblank(2,j) .ne.1) cycle eloop ! skip if incomplete mesh B elem
             qB=mshB%q(1,:,j)

             if(mshA%dx(eid).lt.mshB%dx(j)) then ! A is fine mesh
               xfac = foverlap
             else ! A is coarse mesh
               xfac = 1d0-foverlap
             endif

             if ((x1-y1)*(x1-y2) .le. 0.0) then ! L node of mesh A is inside of mesh B elem
               ! Overlap is between x1 and y2
               ! mshA will remove first half of overlap (from x1 to 0.5*(x1+y2))
               xcut = [x1,x1+xfac*(y2-x1)]
               if(consoverset.eq.1) then 
                 elemInfo(2:3,i) = [xcut(2),x2]
                 mshA%dxcut(i) = x2-xcut(2)
               else
                 elemInfo(2:3,i) = [x1,x2]
               endif
             elseif ((x2-y1)*(x2-y2) .le. 0.0) then ! R node of mesh A is inside of mesh B elem          
               ! Overlap is between y1 and x2
               ! msh A will remove second half of overlap (from 0.5(y1+x2) to x2
               xcut = [x2-xfac*(x2-y1),x2]
               if(consoverset.eq.1) then 
                 elemInfo(2:3,i) = [x1,xcut(1)]
                 mshA%dxcut(i) = xcut(1)-x1
               else
                 elemInfo(2:3,i) = [x1,x2]
               endif
             endif
             lcut = xcut(2)-xcut(1)
             xc = 0.5*(xcut(1)+xcut(2))   ! center of section to be removed
               
             ! Adjust mass matrix 
             if(consoverset.eq.1) then 
               do aa = 1,mshA%ngauss
                 ! get shapefunction from msh A at quad pts of cut section
                 xg = mshA%xgauss(aa)*lcut+xc
                 wtmp = 1d0
                 call shapefunction(mshA%nshp,xg,[x1,x2],wtmp,wtmp,dwtmp)
                 do bb = 1,mshA%nshp
                 do cc = 1,mshA%nshp
                      index1 = (bb-1)*mshA%nshp+cc
                      ! Fix mass matrix
                      mshA%mass(:,index1,eid) = mshA%mass(:,index1,eid) - wtmp(bb)*wtmp(cc)*mshA%wgauss(aa)*lcut
                 enddo ! nshp
  
                 enddo ! nshp
               enddo ! ngauss
  write(*,*) ' '
  write(*,*) 'Mass 2: = ',eid,lcut/mshA%dx(eid),mshA%mass(1,:,eid)
             endif

             cycle iloop
          endif
       enddo eloop
    enddo iloop
  end subroutine fixMassIncompleteElements
  !
end module overset
          
  
