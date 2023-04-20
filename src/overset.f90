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
  subroutine fixFluxIncompleteElements(mshB,mshA,elemInfo,nincomp,consoverset,foverlap,isupg,dt,debug)
    use bases

    ! Subtract half of overlap section from mesh A (stored in elemInfo)
    implicit none
    type(mesh), intent(inout) :: mshA,mshB
    integer, intent(in) :: nincomp,consoverset,isupg,debug
    real*8, intent(inout) :: elemInfo(3,nincomp),foverlap
    real*8, intent(in) :: dt
    !
    integer :: i,j,k,e,nrows,aa,bb,cc,eid,neigh,id,pidA,pidB
    real*8 :: x1,x2,xp1,xp2,f1,f2,y1,y2,yp1,yp2,qA(mshA%nshp),qB(mshB%nshp)
    real*8 :: xrem(2),xcut(2),xc,lcut,xg,vol,flx,qL,qR,fact,xfac
    real*8 :: wtmp(mshA%nshp),dwtmp(mshA%nshp),tmp
    real*8 ::qtmp(mshA%nshp),dqtmp(mshA%nshp),dq,dvol(mshA%nshp),dflx(mshA%nshp)
    real*8 :: qval,dqval,dudt,resid,tau
    !
    ! elemInfo = incomplete elements on mesh A 
    !    (1: local elemID, 2: parent elemID, 3: xcut 1, 4: xcut 2)
    ! msh = mesh info of mesh B
    !
     if(debug.eq.1) then
       write(*,*) ' '
       write(*,*) '--------------------------------------'
       write(*,*) 'Entering fixFluxIncompleteElements'
       write(*,*) '--------------------------------------'
       write(*,*) ' '
       do i=1,mshA%nelem 
        write(*,*) 'A Parents: ',i,mshA%parent(i)
       enddo
       do i=1,mshB%nelem 
        write(*,*) 'B Parents: ',i,mshB%parent(i)
       enddo
     endif
     iloop: do i=1,nincomp       ! Loop through incomplete elem of mesh A
       eid = elemInfo(1,i)
       x1=mshA%xe(1,eid) 
       x2=mshA%xe(2,eid) 
       eloop: do j=1,mshB%nelem  ! Loop through all elem of mesh B
          y1=mshB%xe(1,j)
          y2=mshB%xe(2,j)
          if (x1 > y2 .or. x2 < y1) then ! skip if not overlapping
             cycle eloop
          else
	     if (mshB%iblank(1,j) .ne.1 .and. &
                 mshB%iblank(2,j) .ne.1) cycle eloop ! skip if element blanked
             pidB = mshB%parent(j)
             yp1=mshB%xe(1,pidB)
             yp2=mshB%xe(2,pidB)
             qB=mshB%q(1,:,pidB)

             if(mshA%dx(eid).lt.mshB%dx(j)) then ! A is fine mesh
               xfac = foverlap
             else ! A is coarse mesh
               xfac = 1d0-foverlap
             endif

             if ((x1-y1)*(x1-y2) .le. 0.0) then ! L node of mesh A is inside of mesh B elem
               ! parent element is element on right

               ! Full overlap is between x1 and y2
               ! mshA will remove first section of overlap 
               if(consoverset.eq.1) then 
                 xcut = [x1,x1+xfac*(y2-x1)]
               else
                 xcut = [x1,x1]
                 if(debug.eq.1) write(*,*) 'XCUT = ',xcut,xfac,x1,x2,y1,y2
               endif
               xrem = elemInfo(2:3,i)
               pidA = mshA%parent(eid)
               qA=mshA%q(1,:,pidA) 
               xp1=mshA%xe(1,pidA) 
               xp2=mshA%xe(2,pidA) 
               if(debug.eq.1) then
                 write(*,*) 'L node inside mesh B:'
                 write(*,*) '  Mesh A elem:',eid,x1,x2
                 write(*,*) '  Mesh B elem:',y1,y2
                 write(*,*) '  ElemInfo:',elemInfo(:,i)
                 write(*,*) '  Parent:',pidA,xp1,xp2,qA
                 write(*,*) '  Parent Na, qA= ',wtmp
                 write(*,*) ' ' 
               endif

               ! add intermesh flux from mesh B interior to mesh A L node
               wtmp = 1d0
               call shapefunction(mshA%nshp,xcut(2),[xp1,xp2],wtmp,wtmp,dwtmp)
               call shapefunction(mshB%nshp,xcut(2),[yp1,yp2],qB,qtmp,dqtmp)
               qL = sum(qtmp)
               call shapefunction(mshA%nshp,xcut(2),[xp1,xp2],qA,qtmp,dqtmp)
               qR = sum(qtmp)
               call flux(qL,qR,flx)
               do k = 1,mshA%nshp
                 mshA%rhs(:,k,pidA) = mshA%rhs(:,k,pidA) + wtmp(k)*flx
               enddo
               if(debug.eq.1) write(*,*) 'L Debug fflux 1: ',xcut(2),wtmp,qL,qR,flx

               if(pidA.eq.eid) then ! only do if not using cell agglomeration
                 ! Handle mesh A R node flux
                 wtmp = 1d0
                 call shapefunction(mshA%nshp,x2,[x1,x2],wtmp,wtmp,dwtmp)
                 neigh = mshA%face(2,eid)
                 call shapefunction(mshA%nshp,x2,[xp1,xp2],qA,qtmp,dqtmp)
                 qL = sum(qtmp)
                 call shapefunction(mshA%nshp,x2,mshA%xe(:,neigh),mshA%q(1,:,neigh),qtmp,dqtmp)
                 qR = sum(qtmp)
                 call flux(qL,qR,flx)
                 do k = 1,mshA%nshp
                   mshA%rhs(:,k,pidA) = mshA%rhs(:,k,pidA) - wtmp(k)*flx
                 enddo
                 if(debug.eq.1) write(*,*) 'R Debug fflux 2: ',x2,wtmp,qL,qR,flx
               endif

             elseif ((x2-y1)*(x2-y2) .le. 0.0) then ! R node of mesh A is inside of mesh B elem          
               ! Overlap is between y1 and x2
               ! msh A will remove second half of overlap (from 0.5(y1+x2) to x2
               if(consoverset.eq.1) then 
                 xcut = [x2-xfac*(x2-y1),x2]
                 if(debug.eq.1) write(*,*) 'XCUT = ',xcut,xfac,x1,x2,y1,y2
               else
                 xcut = [x2,x2]
               endif
               xrem = elemInfo(2:3,i)
               pidA = mshA%parent(eid)
               qA=mshA%q(1,:,pidA) 
               xp1=mshA%xe(1,pidA)
               xp2=mshA%xe(2,pidA)
               if(debug.eq.1) then
                 write(*,*) 'R node inside mesh B:'
                 write(*,*) '  Mesh A elem:',eid,x1,x2
                 write(*,*) '  Mesh B elem:',y1,y2
                 write(*,*) '  ElemInfo:',elemInfo(:,i)
                 write(*,*) '  Parent:',pidA,xp1,xp2,qA
                 write(*,*) ' ' 
               endif

               if(pidA.eq.eid) then ! only do if not using cell agglomeration
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
                   mshA%rhs(:,k,pidA) = mshA%rhs(:,k,pidA) + wtmp(k)*flx
                 enddo
                 if(debug.eq.1)  write(*,*) 'L Debug fflux 3: ',x1,wtmp,qL,qR,flx
               endif

               ! add intermesh flux from mesh B interior to mesh A R node
               ! Note we're adding to the parent element pidA, not eid
               wtmp = 1d0
               call shapefunction(mshA%nshp,xcut(1),[xp1,xp2],wtmp,wtmp,dwtmp)
               call shapefunction(mshB%nshp,xcut(1),[y1,y2],qB,qtmp,dqtmp)
               qR = sum(qtmp)
               call shapefunction(mshA%nshp,xcut(1),[xp1,xp2],qA,qtmp,dqtmp)
               qL = sum(qtmp)
               call flux(ql,qr,flx)
               do k = 1,mshA%nshp
                 mshA%rhs(:,k,pidA) = mshA%rhs(:,k,pidA) - wtmp(k)*flx
               enddo
               if(debug.eq.1)  write(*,*) 'R Debug fflux 4: ',xcut(1),wtmp,qL,qR,flx

             endif
             fact = (xrem(2)-xrem(1))/(x2-x1)
               
             dvol = 0d0
             ! Compute volume integral over partial element using parent element
             ! bases functions
             do aa = 1,mshA%ngauss
               ! get shapefunction from msh A at quad pts of remaining element
               xg = mshA%xgauss(aa)*(xrem(2)-xrem(1))+0.5d0*(xrem(2)+xrem(1))
               wtmp = 1d0
               call shapefunction(mshA%nshp,xg,[xp1,xp2],wtmp,wtmp,dwtmp)
               call shapefunction(mshA%nshp,xg,[xp1,xp2],qA,qtmp,dqtmp)
               qval = sum(qtmp)
               dqval = sum(dqtmp)/mshA%dx(pidA)
               call volint(qval,dqval,vol)
               do bb = 1,mshA%nshp
                 ! Volume Integral
                 mshA%rhs(:,bb,pidA) = mshA%rhs(:,bb,pidA) + dwtmp(bb)*vol*(mshA%wgauss(aa)*fact) ! scale gauss weights by length of remaining element parent

                 dvol(bb) = dvol(bb) + dwtmp(bb)*vol*(mshA%wgauss(aa)*fact)
                 if(debug.eq.1) write(*,*) 'Debug volflux:',pidA,bb,xg,dwtmp(bb),vol
                 if(debug.eq.1) write(*,*) '    ',mshA%wgauss(aa),fact,dwtmp(bb)*vol*(mshA%wgauss(aa)*fact)

!                 ! SUPG Terms
!                 if(isupg.eq.1) then
!                   ! get residual
!                   call shapefunction(mshA%nshp,xg,[x1,x2],mshA%qold(1,:,eid),qtmp,dqtmp)
!                   dudt = (qval-sum(qtmp))/dt
!                   if (index(pde_descriptor,'linear_advection') > 0 ) then
!                     resid = dudt + a*dqval ! du/dt + a du/dx
!                     tau = sqrt(a*a/mshA%dxcut(eid)/mshA%dxcut(eid) + 4d0/dt/dt)
!                     tau = a/tau
!                   else if (index(pde_descriptor,'burgers') > 0) then
!                     resid = dudt + qval*dqval ! du/dt + u du/dx
!                     tau = sqrt(qval*qval/mshA%dxcut(eid)/mshA%dxcut(eid) + 4d0/dt/dt)
!                     tau = qval/tau
!                   endif
!                   if(qval.gt.1.5d0) then 
!                   write(*,*) ' '
!                   write(*,*) ' eid = ',eid
!                   write(*,*) '   u,uold,dt = ',qval, sum(qtmp),dt
!                   write(*,*) '   dudt, dqval =',dudt, dqval
!                   write(*,*) '   resid = ',resid
!                   write(*,*) '   tau = ',tau
!                   write(*,*) '   supg = ',tau*resid
!                   endif
!                   mshA%rhs(:,bb,eid) = mshA%rhs(:,bb,eid) - dwtmp(bb)*tau*resid*(mshA%wgauss(aa)*fact)
!                 endif ! supg
               enddo ! nshp
             enddo ! ngauss
             if(debug.eq.1) write(*,*) '  Total cut volflux: ',dvol,mshA%rhs(:,:,pidA)
             cycle iloop
          endif
       enddo eloop
    enddo iloop
     if(debug.eq.1) then
       write(*,*) ' '
       write(*,*) '--------------------------------------'
       write(*,*) 'Exiting fixFluxIncompleteElements'
       write(*,*) '--------------------------------------'
       write(*,*) ' '
     endif
  end subroutine fixFluxIncompleteElements
  !
  subroutine fixMassIncompleteElements(mshB,mshA,elemInfo,nincomp,consoverset,foverlap,debug)
    use bases

    ! Subtract half of overlap section from mesh A (stored in elemInfo)
    implicit none
    type(mesh), intent(inout) :: mshA,mshB
    integer, intent(in) :: nincomp,consoverset,debug
    real*8, intent(inout) :: elemInfo(3,nincomp),foverlap
    !
    integer :: i,j,k,e,nrows,aa,bb,cc,eid,index1,pidA
    real*8 :: x1,x2,f1,f2,y1,y2,qA(mshA%nshp),qB(mshB%nshp),xp1,xp2
    real*8 :: xcut(2),xc,lcut,xg,xfac
    real*8 :: wtmp(mshA%nshp),dwtmp(mshA%nshp)
    !
    ! elemInfo = incomplete elements on mesh A
    ! msh = mesh info of mesh B
    !
     if(debug.eq.1) then
       write(*,*) ' '
       write(*,*) '--------------------------------------'
       write(*,*) 'Entering fixMassIncompleteElements'
       write(*,*) '--------------------------------------'
       write(*,*) ' '
     endif
     iloop: do i=1,nincomp       ! Loop through incomplete elem of mesh A
       eid = elemInfo(1,i)
       x1=mshA%xe(1,eid) 
       x2=mshA%xe(2,eid) 
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
                 if(debug.eq.1) write(*,*) 'Elem ',i,'cut ratio = ',mshA%dxcut(i)/mshA%dx(i)
               else
                 elemInfo(2:3,i) = [x1,x2]
               endif
               if((mshA%dxcut(i)/mshA%dx(i).le.0.1d0).and.(consoverset.eq.1)) then 
!               if((consoverset.eq.1)) then 
                 pidA = mshA%face(2,eid)
                 mshA%child(pidA) = eid
                 if(debug.eq.1) write(*,*) 'MERGING CELL ',eid,' AND ',pidA
               else
                 pidA = eid
               endif
               mshA%parent(eid) = pidA
               xp1=mshA%xe(1,pidA)
               xp2=mshA%xe(2,pidA)
               if(debug.eq.1) then
                 write(*,*) 'L node inside mesh B:'
                 write(*,*) '  Mesh A elem:',eid,x1,x2
                 write(*,*) '  Mesh B elem:',y1,y2
                 write(*,*) '  xcut:',xcut(1),xcut(2)
                 write(*,*) '  ElemInfo:',elemInfo(:,i)
                 write(*,*) '  Parent:',pidA,xp1,xp2
                 write(*,*) '  Parents Child:',mshA%child(pidA)
                 write(*,*) ' ' 
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
               if((mshA%dxcut(i)/mshA%dx(i).le.0.1d0).and.(consoverset.eq.1)) then 
!               if((consoverset.eq.1)) then 
                 pidA = mshA%face(1,eid)
                 mshA%child(pidA) = eid
                 if(debug.eq.1) write(*,*) 'MERGING CELL ',eid,' AND ',pidA
               else
                 pidA = eid
               endif
               mshA%parent(eid) = pidA
               xp1=mshA%xe(1,pidA)
               xp2=mshA%xe(2,pidA)
               if(debug.eq.1) then
                 write(*,*) 'R node inside mesh B:'
                 write(*,*) '  Mesh A elem:',eid,x1,x2
                 write(*,*) '  Mesh B elem:',y1,y2
                 write(*,*) '  xcut:',xcut(1),xcut(2)
                 write(*,*) '  ElemInfo:',elemInfo(:,i)
                 write(*,*) '  Parent:',pidA,xp1,xp2
                 write(*,*) '  Parents Child:',mshA%child(pidA)
                 write(*,*) ' ' 
               endif
             endif ! L or R side
             if(debug.eq.1) then
               write(*,*) '  Original Mass:',mshA%mass(1,:,pidA)
               write(*,*) ' ' 
             endif
               
             ! Adjust mass matrix 
             if(consoverset.eq.1) then 
               if(pidA.ne.eid) then ! if agglomerated cell
                 ! add additional mass over full child cell
                 do aa = 1,mshA%ngauss
                   xg = mshA%xgauss(aa)*mshA%dx(eid)+0.5d0*(mshA%xe(1,eid)+mshA%xe(2,eid))
                   wtmp = 1d0
                   call shapefunction(mshA%nshp,xg,[xp1,xp2],wtmp,wtmp,dwtmp)
                   do bb = 1,mshA%nshp
                     do cc = 1,mshA%nshp
                       index1 = (bb-1)*mshA%nshp+cc
                       ! Fix mass matrix
                       mshA%mass(:,index1,pidA) = mshA%mass(:,index1,pidA) + wtmp(bb)*wtmp(cc)*mshA%wgauss(aa)*mshA%dx(eid)
                     enddo ! nshp
                   enddo ! nshp
                 enddo ! ngauss
               endif

               ! subtract child cell's cut portion
               lcut = xcut(2)-xcut(1)
               xc = 0.5*(xcut(1)+xcut(2))   ! center of section to be removed
               do aa = 1,mshA%ngauss
                 ! get shapefunction from msh A at quad pts of cut section
                 xg = mshA%xgauss(aa)*lcut+xc
                 wtmp = 1d0
                 call shapefunction(mshA%nshp,xg,[xp1,xp2],wtmp,wtmp,dwtmp)
                 do bb = 1,mshA%nshp
                 do cc = 1,mshA%nshp
                      index1 = (bb-1)*mshA%nshp+cc
                      ! Fix mass matrix
                      mshA%mass(:,index1,pidA) = mshA%mass(:,index1,pidA) - wtmp(bb)*wtmp(cc)*mshA%wgauss(aa)*lcut
                 enddo ! nshp
  
                 enddo ! nshp
               enddo ! ngauss
             endif
             if(debug.eq.1) then
               write(*,*) '  Modified Mass:',mshA%mass(1,:,pidA)
               write(*,*) ' ' 
             endif

             cycle iloop
          endif
       enddo eloop
    enddo iloop
     if(debug.eq.1) then
       write(*,*) ' '
       write(*,*) '--------------------------------------'
       write(*,*) 'Exiting fixMassIncompleteElements'
       write(*,*) '--------------------------------------'
       write(*,*) ' '
     endif
  end subroutine fixMassIncompleteElements
  !
  subroutine projectChild(msh,elemInfo,nincomp,vec)
    use bases

    implicit none
    type(mesh),intent(inout) :: msh
    integer :: i,j,k,eid,pid
    integer,intent(in) :: nincomp
    real*8, intent(inout) :: vec(msh%nfields,msh%nshp,msh%nelem)
    real*8,dimension(msh%nshp*msh%nshp) :: L,U
    real*8, intent(in) :: elemInfo(3,nincomp)
    real*8,dimension(msh%nshp) :: y
    real*8 :: x1,x2,xloc,qvals(msh%nshp),dqvals(msh%nshp),xp1,xp2
    real*8 :: f(msh%nshp),dx,tmp

    do i = 1,nincomp
      eid = elemInfo(1,i)
      pid = msh%parent(eid)
      xp1=msh%xe(1,pid)
      xp2=msh%xe(2,pid)

      if(pid.ne.eid) then
        if(shptype.eq.'legendre') then 
          ! Use projection to set the IC
          ! int(NaNbub) = int(Na q_extrap)
          ! ub = M^-1 int(Na q_extrap)
          dx = msh%dx(eid)
          f = 0.0d0
          do j = 1,msh%ngauss
            ! extrapolate q value based on parent bases
            xloc = (msh%xgauss(j)+0.5d0)*dx + msh%xe(1,eid)
            call shapefunction(msh%nshp,xloc,[xp1,xp2],vec(1,:,pid),qvals,dqvals)
            tmp = SUM(qvals)
            qvals = 1d0
            call shapefunction(msh%nshp,msh%xgauss(j),[-0.5d0,0.5d0],qvals,qvals,dqvals)
            do k = 1,msh%nshp
              f(k) = f(k) + msh%wgauss(j)*msh%dx(pid)*tmp*qvals(k)
            enddo
          enddo

          ! solve the elementary system for u
          call lu(msh%mass(1,:,pid),msh%nshp,L,U)
          call forwprop(L,f,msh%nshp,y)
          call backprop(U,y,msh%nshp,vec(1,:,eid))
        else if(shptype.eq.'lagrange') then        
          do j = 1,msh%nshp
            ! compute q value at nodal point using parent element Q values
            xloc = msh%x(msh%e2n(j,eid))
            call shapefunction(msh%nshp,xloc,[xp1,xp2],vec(1,:,pid),qvals,dqvals)
            vec(1,j,eid) = SUM(qvals)
          enddo ! loop over shape functions
        endif ! shp type
      endif ! if merged cell
    enddo ! loop over nincomp
  end subroutine projectChild
  !
end module overset
          
  
