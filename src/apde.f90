module pde
  use code_types
  use bases
  character*20 :: pde_descriptor
  integer, save :: pde_set=0
  real*8 :: a,mu,qin,qout,porder
contains
  subroutine set_type(pdetype,wavespeed,muinf,qA,qB,p)
    implicit none
    character*(*) :: pdetype
    real*8, optional, intent(in) :: wavespeed,muinf,qA,qB
    integer :: p
    write(pde_descriptor,'(A20)') pdetype
    pde_descriptor=trim(adjustl(pde_descriptor))
    if (pde_descriptor .ne. 'linear_advection' .and. & 
        pde_descriptor .ne. 'burgers' ) then
       write(6,*) 'Only linear advection and burgers are implemented'
       stop
    endif
    pde_set=1
    a=wavespeed
    mu=muinf
    qin = qA
    qout = qB
    porder=p
  end subroutine set_type
  !
  subroutine initqleg(msh,q,t)
    use code_types
    use bases
    type(mesh), intent(inout)::msh
    integer i,j,k
    real*8, intent(inout) :: q(msh%nfields,msh%nshp,msh%nelem)
    real*8,dimension(msh%nshp*msh%nshp) :: L,U
    real*8,dimension(msh%nshp) :: y
    real*8 :: qvals(msh%nshp),dqvals(msh%nshp),f(msh%nshp),tmp,xloc,dx,t
    !
    ! Use projection to set the IC
    ! int(NaNbub) = int(Na y(x,0))
    ! ub = M^-1 int(Na y(x,0))
    do i = 1,msh%nelem
      dx = msh%xe(2,i)-msh%xe(1,i)
      f = 0.0d0
      do j = 1,msh%ngauss
        xloc = (msh%xgauss(j)+0.5d0)*dx + msh%xe(1,i)
        if(index(pde_descriptor,'burgers') > 0) then ! burgers
          tmp = 1d0-tanh((xloc+0.5-t)/(2*mu))
        else ! A-D
          xloc = mod(1d0+xloc-t*a,2d0)-1d0
          tmp = 0.1d0*exp(-20d0*xloc*xloc)
          if (abs(xloc) > 0.8d0) tmp=0d0
        endif
        qvals = 1d0
        call shapefunction(msh%nshp,msh%xgauss(j),[-0.5d0,0.5d0],qvals,qvals,dqvals)
        do k = 1,msh%nshp
          f(k) = f(k) + msh%wgauss(j)*msh%dx(i)*tmp*qvals(k)
        enddo
      enddo

      ! solve the elementary system for u
      call lu(msh%mass(1,:,i),msh%nshp,L,U)
      call forwprop(L,f,msh%nshp,y)
      call backprop(U,y,msh%nshp,q(1,:,i))

    enddo
  end subroutine initqleg
  !
  subroutine initq(x,q,t)
    implicit none
    real*8, intent(in) :: x,t
    real*8, intent(out) :: q(:)
    real*8 :: xx

    if ((index(pde_descriptor,'linear_advection') .le. 0) .and. &
        (index(pde_descriptor, 'burgers') .le. 0)) then
       write(6,*) 'problem type not implemented'
       call exit(1)
    else     
       if(index(pde_descriptor,'burgers') > 0) then ! burgers
         q=1d0-tanh((x+0.5-t)/(2*mu))
       else ! A-D
         xx = mod(1d0+x-t*a,2d0)-1d0 ! assume periodic domain size 2
         q=exp(-20*xx*xx)
         if (abs(x) > 0.8d0) q=0d0
       endif

!       q = 1d0
!       if (x < 1d0) q=0d0
!       if (x > 5d0) q=0d0
!       q=exp(-0.01*x*x)
!       if (abs(x) > 40d0) q=0d0
    endif
  end subroutine initq
  !
  ! New flux routine
  ! Outputs entire bounday integral vector
  subroutine flux2(mshA,eidA,mshB,eidB,xloc,flx)
    !
    implicit none
    !
    type(mesh), intent(inout) :: mshA, mshB
    real*8, intent(in) :: xloc
    real*8, intent(out) :: flx(mshA%nshp)
    integer :: eidA, eidB
    real*8 :: qAvals(mshA%nshp),dqAvals(mshA%nshp),wAvals(mshA%nshp),dwAvals(mshA%nshp)
    real*8 :: qA,dqA,wA,dwA
    real*8 :: qBvals(mshA%nshp),dqBvals(mshA%nshp),wBvals(mshA%nshp),dwBvals(mshA%nshp)
    real*8 :: qB,dqB,wB,dwB
    real*8 :: C11,C12,flxa(mshA%nshp),flxd(mshA%nshp),nA,nB
    !
    flx = 0d0
    !
    ! Get Interior State (A)
    wAvals = 1d0
    call shapefunction(mshA%nshp,xloc,mshA%xe(:,eidA),wAvals,wAvals,dwAvals)
    wA = sum(wAvals)
    dwA = sum(dwAvals)/mshA%dx(eidA)
    call shapefunction(mshA%nshp,xloc,mshA%xe(:,eidA),mshA%q(1,:,eidA),qAvals,dqAvals)
    qA = sum(qAvals)
    dqA = sum(dqAvals)/mshA%dx(eidA)
    if((mshA%iBC(eidA).eq.1).and.(xloc.eq.mshA%xe(1,eidA))) then ! inflow
      qA = qin
      dqA = 0d0
    else if((mshA%iBC(eidA).eq.-1).and.(xloc.eq.mshA%xe(2,eidA))) then ! outflow
      qA = qout
      dqA = 0d0
    endif

    ! Get Exterior State (B)
    wAvals = 1d0
    call shapefunction(mshB%nshp,xloc,mshB%xe(:,eidB),wBvals,wBvals,dwBvals)
    wB = sum(wBvals)
    dwB = sum(dwBvals)/mshB%dx(eidB)
    call shapefunction(mshB%nshp,xloc,mshB%xe(:,eidB),mshB%q(1,:,eidB),qBvals,dqBvals)
    qB = sum(qBvals)
    dqB = sum(dqBvals)/mshB%dx(eidB)
    if((mshB%iBC(eidB).eq.1).and.(xloc.eq.mshB%xe(1,eidB))) then ! inflow
      qB = qin
      dqB = 0d0
    else if((mshB%iBC(eidB).eq.-1).and.(xloc.eq.mshB%xe(2,eidB))) then ! outflow
      qB = qout
      dqB = 0d0
    endif

    ! computing normals
    ! ie figuring out if cell A is L or R cell
    if(sum(mshA%xe(:,eidA)).lt.sum(mshB%xe(:,eidB))) then
      nA = 1d0
    else
      nA = -1d0
    endif
    nB = -nA

    ! Advective fluxes
    if (index(pde_descriptor,'linear_advection') > 0 ) then
      ! LDG method
      C11 = 0.0d0
      C12 = 0.5d0 ! 0.5 dot n-
      flxa=wAvals*(a*(0.5d0*(qA+qB) + C12*(qA-qB)))
    else if (index(pde_descriptor,'burgers') > 0) then
      ! LDG method
      C11 = 0.0d0
      C12 = 0.5d0 ! 0.5 dot n-
      flxa = wAvals*(0.5d0*(0.5d0*(qA*qA+qB*qB)-0.5d0*(qA+qB)*(qB-qA)))

      ! IP or BR2 symmetric adv flux
      ! flxa = wAvals*0.5d0*(0.5d0*(qA*qA+qB*qB))
    endif
   
    ! Diffusive Fluxes
    if (index(pde_descriptor,'linear_advection') > 0 ) then
      ! LDG method
      C11 = 0.0d0
      C12 = 0.5d0 ! 0.5 dot n-

      flxd = wAvals*mu*(0.5d0*(dqA+dqB) + C11*(qA-qB) - C12*(dqA-dqB)) ! LDG method
    else if (index(pde_descriptor,'burgers') > 0) then
      ! Advective Flux
      ! Lax Friedrich
      C11 = 0.0d0
      C12 = 0.5d0 ! 0.5 dot n-
      flxa = wAvals*0.5d0*(0.5d0*(qA*qA+qB*qB)-0.5d0*(qA+qB)*(qB-qA))

      ! Arnold/Shabazi Interior Penalty:
      C11 = (porder+1d0)*(porder+1d0)/2d0 ! Shabazi penalty constant
      C11 = C11/(mshA%dx(eidA))
       
      flxd = 0.5d0*(dqAvals+dqBvals)*(wA*nA+wB*nB)
      flxd = flxd + 0.5d0*(dwAvals+dwBvals)*(qA*nA+qB*nB)
      flxd = flxd - wAvals*C11*(qA*nA+qB*nB)*(wA*nA+wB*nB) ! penalty term

      ! LDG according to Persson DG School
      ! works best so far except on fine meshes & p>4
!      flxd = wAvals*mu*(0.5d0*(dql+dqr) + C11*(ql-qr) - C12*(dql-dqr)) ! LDG method

      ! Westhaven and Warburton:

      ! Symmetric flux
!      flxd = wAvals*0.5d0*mu*(dqA+dqB)
    endif
 
    ! Final Flux
    flx=flxa-flxd
  end subroutine flux2

  subroutine flux(ql,qr,dql,dqr,flx)
    !
    implicit none
    !
    real*8, intent(in) :: ql,dql
    real*8, intent(in) :: qr,dqr
    real*8, intent(out) :: flx
    real*8 :: flxa,flxd,vfl,vfr,penalty,tm
    real*8 :: C11,C12
    !
    if (index(pde_descriptor,'linear_advection') > 0 ) then
      ! Lax Friedrichs Flux
      ! flxa=0.5d0*(a*(ql+qr)-abs(a)*(qr-ql))

      C11 = 0.0d0
      C12 = 0.5d0 ! 0.5 dot n-

      flxa=a*(0.5d0*(ql+qr) + C12*(ql-qr) )
      flxd = mu*(0.5d0*(dql+dqr) + C11*(ql-qr) - C12*(dql-dqr)) ! LDG method

      flx=flxa-flxd
    else if (index(pde_descriptor,'burgers') > 0) then
      ! Advective Flux
      ! Lax Friedrich
      C11 = 0.0d0
      C12 = 0.5d0 ! 0.5 dot n-
      flxa = 0.5d0*(0.5d0*(ql*ql+qr*qr)-0.5d0*(ql+qr)*(qr-ql))

      ! LDG according to Persson DG School
      ! works best so far except on fine meshes & p>4
      flxd = mu*(0.5d0*(dql+dqr) + C11*(ql-qr) - C12*(dql-dqr)) ! LDG method
      ! Westhaven and Warburton:
      
      ! Penalty from Brazell matlab code
!      tm = 5 ! p+1, testing for p=5
!      penalty = tm/0.0078125d0 ! denom is dx/2
!      vfl = mu*(qL-qR)/2d0
!      vfr = mu*(qR-qL)/2d0
!      flxd = mu*(0.5d0*(dql+dqr)) -penalty*(vfl-vfr)
      ! Symmetric flux
!      flxd = 0.5d0*mu*(dql+dqr)
      flx=flxa-flxd
    endif
  end subroutine flux
  !
  subroutine volint(qtmp,dqtmp,vol)
    !
    implicit none
    !
    real*8, intent(in) :: qtmp,dqtmp
    real*8, intent(out) :: vol
    !
    if (index(pde_descriptor,'linear_advection') > 0 ) then
     vol=a*qtmp - mu*dqtmp
    else if (index(pde_descriptor,'burgers') > 0) then
     vol=0.5*qtmp*qtmp - mu*dqtmp
    endif
  end subroutine volint
  ! 
end module pde
