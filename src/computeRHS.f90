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
  !
  do i = 1,msh%nelem
    e1=msh%face(1,i)
    e2=msh%face(2,i)
    ! Calculate the volume integrals
    ! For Fringe elements, do the volume integral but not fluxes
    do j = 1,msh%ngauss
      call shapefunction(msh%nshp,msh%xgauss(j),[-0.5d0,0.5d0],msh%q(1,:,i),qvals,dqvals)
      qtmp = SUM(qvals)
      do k = 1,msh%nshp
         call volint(qtmp,vol)
         vol = vol*msh%dshp(j,k)*msh%wgauss(j)
         msh%rhs(1,k,i) = msh%rhs(1,k,i) + vol 
      enddo 
    enddo
    !
    ! Calculate the fluxes within current mesh
    ! For fringe elements where flux comes from other mesh, 
    ! intermesh fluxes done in overset routines
    if (e1.ne.e2) then

        ! Left flux boundary        
        call shapefunction(msh%nshp,-0.5d0,[-0.5d0,0.5d0],[1d0,1d0],qvals,dqvals)
        w = qvals
        call shapefunction(msh%nshp,msh%xe(2,e1),msh%xe(:,e1),msh%q(1,:,e1),qvals,dqvals)
        ql = sum(qvals)
        call shapefunction(msh%nshp,msh%xe(1,i),msh%xe(:,i),msh%q(1,:,i),qvals,dqvals)
        qr = sum(qvals)
        call flux(ql,qr,flx)
        if ((msh%iblank(msh%e2n(2,e1)) > 0).and.(e1.ne.i)) then 
          do j = 1,msh%nshp
            msh%rhs(:,j,i) = msh%rhs(:,j,i) + w(j)*flx
          enddo
        endif
        
        ! Right flux boundary
        call shapefunction(msh%nshp,0.5d0,[-0.5d0,0.5d0],[1d0,1d0],qvals,dqvals)
        w = qvals
        call shapefunction(msh%nshp,msh%xe(2,i),msh%xe(:,i),msh%q(1,:,i),qvals,dqvals)
        ql = sum(qvals)
        call shapefunction(msh%nshp,msh%xe(1,e2),msh%xe(:,e2),msh%q(1,:,e2),qvals,dqvals)
        qr = sum(qvals)
        call flux(ql,qr,flx)
        if ((msh%iblank(msh%e2n(1,e2)) > 0).and.(e2.ne.i)) then 
          do j = 1,msh%nshp
            msh%rhs(:,j,i) = msh%rhs(:,j,i) - w(j)*flx
          enddo
        endif
    endif
  enddo
  !
  !iout=iout+1
  !do i=1,msh%nelem
  ! write(30+iout,*) msh%rhs(1,1,i),msh%iblank(msh%e2n(1,i)),&
  !                  msh%iblank(msh%e2n(2,i))
  !enddo
end subroutine computeRHS
     
     
