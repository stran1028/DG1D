subroutine computeRHS(msh)
  use code_types
  use pde
  use bases
  implicit none
  !
  type(mesh), intent(inout) :: msh
  integer :: i,j,k,e1,e2
  real*8 :: ql,qr,flx,qtmp,qvals(msh%nshp),dqvals(msh%nshp),vol,w(msh%nshp)
  integer, save :: iout=0
  !
  ! this is p0 implementation for
  ! single field
  !
  ! rhs = (msh%nfields,msh%nshp,msh%nelem)
  msh%rhs=0d0
  msh%rhsV=0d0
  msh%rhsF=0d0
  !
  do i = 1,msh%nelem
    e1=msh%face(1,i)
    e2=msh%face(2,i)

    ! Only look at interior mesh elements, save fringes for later
    if(minval(msh%iblank(msh%e2n(:,i))).eq.1d0) then 
      ! Calculate the volume integrals
      do j = 1,msh%ngauss
        call shapefunction(msh%nshp,msh%xgauss(j),[-0.5d0,0.5d0],msh%q(1,:,i),qvals,dqvals)
        qtmp = SUM(qvals)
        do k = 1,msh%nshp
           call volint(qtmp,vol)
           vol = vol*msh%dshp(j,k)*msh%wgauss(j)
           msh%rhs(1,k,i) = msh%rhs(1,k,i) + vol 
           msh%rhsV(1,k,i) = msh%rhsV(1,k,i) + vol 
        enddo 
      enddo
      !
      ! Calculate the fluxes within current mesh
      if (e1.ne.e2) then
        ! Left flux boundary        
        call shapefunction(msh%nshp,-0.5d0,[-0.5d0,0.5d0],[1d0,1d0],qvals,dqvals)
        w = qvals
        call shapefunction(msh%nshp,msh%xe(2,e1),msh%xe(:,e1),msh%q(1,:,e1),qvals,dqvals)
        ql = sum(qvals)
        call shapefunction(msh%nshp,msh%xe(1,i),msh%xe(:,i),msh%q(1,:,i),qvals,dqvals)
        qr = sum(qvals)
        call flux(ql,qr,flx)
        do j = 1,msh%nshp
            msh%rhs(:,j,i) = msh%rhs(:,j,i) + w(j)*flx
            msh%rhsF(:,j,i) = msh%rhsF(:,j,i) + w(j)*flx
        enddo
        
        ! Right flux boundary
        call shapefunction(msh%nshp,0.5d0,[-0.5d0,0.5d0],[1d0,1d0],qvals,dqvals)
        w = qvals
        call shapefunction(msh%nshp,msh%xe(2,i),msh%xe(:,i),msh%q(1,:,i),qvals,dqvals)
        ql = sum(qvals)
        call shapefunction(msh%nshp,msh%xe(1,e2),msh%xe(:,e2),msh%q(1,:,e2),qvals,dqvals)
        qr = sum(qvals)
        call flux(ql,qr,flx)
        do j = 1,msh%nshp
            msh%rhs(:,j,i) = msh%rhs(:,j,i) - w(j)*flx
            msh%rhsF(:,j,i) = msh%rhsF(:,j,i) - w(j)*flx
        enddo
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
     
     
