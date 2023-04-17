subroutine computeRHS(msh,isupg,dt)
  use code_types
  use pde
  use bases
  implicit none
  !
  type(mesh), intent(inout) :: msh
  real*8, intent(in) :: dt
  integer :: i,j,k,e1,e2,isupg,cid
  real*8 :: ql,qr,flx,qtmp,dqtmp,qvals(msh%nshp),dqvals(msh%nshp),vol,w(msh%nshp)
  real*8 :: dudt,resid,tau
  integer, save :: iout=0
  !
  ! this is p0 implementation for
  ! single field
  !
  ! rhs = (msh%nfields,msh%nshp,msh%nelem)
  msh%rhs=0d0
  !
  do i = 1,msh%nelem
    e1=msh%face(1,i)
    e2=msh%face(2,i)

    ! Only look at interior mesh elements, save fringes for later
    if(minval(msh%iblank(:,i)).eq.1) then 
      ! Calculate the volume integrals
      do j = 1,msh%ngauss
        call shapefunction(msh%nshp,msh%xgauss(j),[-0.5d0,0.5d0],msh%q(1,:,i),qvals,dqvals)
        qtmp = SUM(qvals)
        dqtmp = SUM(dqvals)/msh%dx(i)
        do k = 1,msh%nshp
           call volint(qtmp,vol)
           vol = vol*msh%dshp(j,k)*msh%wgauss(j)
           msh%rhs(1,k,i) = msh%rhs(1,k,i) + vol 

!           if(isupg.eq.1) then
!             ! get residual
!             call shapefunction(msh%nshp,msh%xgauss(j),[-0.5d0,0.5d0],msh%qold(1,:,i),qvals,dqvals)
!             dudt = (qtmp-sum(qvals))/dt
!             if (index(pde_descriptor,'linear_advection') > 0 ) then
!               resid = dudt + a*dqtmp ! du/dt + a du/dx
!               tau = sqrt(a*a/msh%dxcut(i)/msh%dxcut(i) + 4d0/dt/dt)
!               tau = a/tau
!             else if (index(pde_descriptor,'burgers') > 0) then
!               resid = dudt + qtmp*dqtmp ! du/dt + u du/dx
!               tau = sqrt(qtmp*qtmp/msh%dxcut(i)/msh%dxcut(i) + 4d0/dt/dt)
!               tau = qtmp/tau
!              endif
!!               if(i.eq.10) then 
!!                       write(*,*) 'wgauss(j) = ',msh%wgauss(j)
!!                       write(*,*) 'dudt,dudx = ',dudt,dqtmp
!!                       write(*,*) 'resid,tau = ',resid,tau
!!             endif
!             msh%rhs(1,k,i) = msh%rhs(1,k,i) - tau*resid*msh%dshp(j,k)*msh%wgauss(j)
!           endif ! supg
        enddo ! shape
      enddo ! gauss
      !
      ! Calculate the fluxes within current mesh
      cid = msh%child(i)
      if (e1.ne.e2) then
        ! Left flux boundary       
        if(cid.ge.i) then          
          qvals = 1d0
          call shapefunction(msh%nshp,-0.5d0,[-0.5d0,0.5d0],qvals,w,dqvals)
          call shapefunction(msh%nshp,msh%xe(2,e1),msh%xe(:,e1),msh%q(1,:,e1),qvals,dqvals)
          ql = sum(qvals)
          call shapefunction(msh%nshp,msh%xe(1,i),msh%xe(:,i),msh%q(1,:,i),qvals,dqvals)
          qr = sum(qvals)
          call flux(ql,qr,flx)
          do j = 1,msh%nshp
              msh%rhs(:,j,i) = msh%rhs(:,j,i) + w(j)*flx
          enddo
        endif
        
        if(cid.le.i) then 
          ! Right flux boundary
          qvals = 1d0
          call shapefunction(msh%nshp,0.5d0,[-0.5d0,0.5d0],qvals,w,dqvals)
          call shapefunction(msh%nshp,msh%xe(2,i),msh%xe(:,i),msh%q(1,:,i),qvals,dqvals)
          ql = sum(qvals)
          call shapefunction(msh%nshp,msh%xe(1,e2),msh%xe(:,e2),msh%q(1,:,e2),qvals,dqvals)
          qr = sum(qvals)
          call flux(ql,qr,flx)
          do j = 1,msh%nshp
              msh%rhs(:,j,i) = msh%rhs(:,j,i) - w(j)*flx
          enddo
        endif
      endif ! e1 ne e2
    endif ! iblank 

  enddo
  !
  !iout=iout+1
  !do i=1,msh%nelem
  ! write(30+iout,*) msh%rhs(1,1,i),msh%iblank(msh%e2n(1,i)),&
  !                  msh%iblank(msh%e2n(2,i))
  !enddo
end subroutine computeRHS
     
     
