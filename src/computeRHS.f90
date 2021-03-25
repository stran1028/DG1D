subroutine computeRHS(msh)
  use code_types
  use pde
  use bases
  implicit none
  !
  type(mesh), intent(inout) :: msh
  integer :: i,j,k,e1,e2
  real*8 :: ql,qr,flx,qtmp,qvals(msh%nshp),dqvals(msh%nshp),vol
  integer, save :: iout=0
  !
  ! this is p0 implementation for
  ! single field
  !
  ! rhs = (msh%nfields,msh%nshp,msh%nelem)
  msh%rhs=0d0
  msh%rhsv=0d0
  msh%rhsf=0d0
  !
  do i = 1,msh%nelem
    ! Calculate the volume integrals
    do j = 1,msh%ngauss
      call shapefunction(msh%nshp,msh%xgauss(j),[-0.5d0,0.5d0],msh%q(1,:,i),qvals,dqvals)
      qtmp = SUM(qvals)
      do k = 1,msh%nshp
         call volint(qtmp,vol)
         vol = vol*msh%dshp(j,k)*msh%detJ(i)*msh%wgauss(j)
         msh%rhs(1,k,i) = msh%rhs(1,k,i) + vol 
         msh%rhsv(1,k,i) = msh%rhsv(1,k,i) + vol 
      enddo 
!      if(abs(i-50).lt.10) write(*,*) 'qtmp = ',i,qtmp,vol,msh%rhsv(1,:,i)
    enddo
    !
    ! Calculate the fluxes
    e1=msh%face(1,i)
    e2=msh%face(2,i)
    !write(*,*) 'element i,e1,e2: ',i,e1,e2
    if (e1 .ne. e2) then
        !if (e1==1 .or. e2 == 1) write(6,*) 'i=',i
        !msh%rhs(:,1,e1)=msh%rhs(:,1,e1)+flx
        !msh%rhs(:,1,e2)=msh%rhs(:,1,e2)-flx

        ! Left flux boundary
        call shapefunction(msh%nshp,msh%xe(2,e1),msh%xe(:,e1),msh%q(1,:,e1),qvals,dqvals)
        ql = sum(qvals)
        call shapefunction(msh%nshp,msh%xe(1,i),msh%xe(:,i),msh%q(1,:,i),qvals,dqvals)
        qr = sum(qvals)
        call flux(ql,qr,flx)
        !write(*,*) '  l flux: ',ql,qr,flx
        if (msh%iblank(msh%e2n(2,e1)) > 0) msh%rhs(:,1,i)=msh%rhs(:,1,i)+flx*msh%detJ(i)
        if (msh%iblank(msh%e2n(2,e1)) > 0) msh%rhsf(:,1,i)=msh%rhsf(:,1,i)+flx*msh%detJ(i)
        
        ! Right flux boundary
!        ql = msh%q(1,msh%nshp,i)
!        qr = msh%q(1,1,e2)
        call shapefunction(msh%nshp,msh%xe(2,i),msh%xe(:,i),msh%q(1,:,i),qvals,dqvals)
        ql = sum(qvals)
        call shapefunction(msh%nshp,msh%xe(1,e2),msh%xe(:,e2),msh%q(1,:,e2),qvals,dqvals)
        qr = sum(qvals)
        call flux(ql,qr,flx)
        !write(*,*) '  r flux: ',ql,qr,flx
        if (msh%iblank(msh%e2n(1,e2)) > 0) msh%rhs(:,msh%nshp,i)=msh%rhs(:,msh%nshp,i)-flx*msh%detJ(i)
        if (msh%iblank(msh%e2n(1,e2)) > 0) msh%rhsf(:,msh%nshp,i)=msh%rhsf(:,msh%nshp,i)-flx*msh%detJ(i)
    endif
    !! Need to do anything special for periodic? XXX
  enddo
  !
  !iout=iout+1
  !do i=1,msh%nelem
  ! write(30+iout,*) msh%rhs(1,1,i),msh%iblank(msh%e2n(1,i)),&
  !                  msh%iblank(msh%e2n(2,i))
  !enddo
end subroutine computeRHS
     
     
