subroutine computeRHS(msh)
  use code_types
  use pde
  implicit none
  !
  type(mesh), intent(inout) :: msh
  integer :: i,j,k,e1,e2
  real*8 :: ql,qr,flx,qtmp,qvals(2),dqvals(2),vol
  integer, save :: iout=0
  !
  ! this is p0 implementation for
  ! single field
  !
  ! rhs = (msh%nfields,msh%nshp,msh%nelem)
  msh%rhs=0d0
  !
  do i = 1,msh%nelem

    ! Calculate the volume integrals
    do j = 1,msh%ngauss
      call shapefunction(msh%nshp,msh%xgauss(j),[-0.5d0,0.5d0],msh%q(1,:,i),qvals,dqvals)
      qtmp = SUM(qvals)
      do k = 1,msh%nshp
         vol = 0d0
         call volint(qtmp,vol)
         vol = vol*msh%dshp(j,k)*msh%detJ(i)*msh%wgauss(j)
         msh%rhs(1,k,i) = msh%rhs(1,k,i) + vol 
      enddo 
    enddo
    if(i.eq.50)  write(*,*) '  q1,q2 = ',msh%q(1,1,i),msh%q(1,msh%nshp,i)
    if(i.eq.50)  write(*,*) '  volume integral:', msh%rhs(1,1,i),msh%rhs(1,2,i)
    !
    ! Calculate the fluxes
    e1=msh%face(1,i)
    e2=msh%face(2,i)
    if (e1 .ne. e2) then
        !if (e1==1 .or. e2 == 1) write(6,*) 'i=',i
        !msh%rhs(:,1,e1)=msh%rhs(:,1,e1)+flx
        !msh%rhs(:,1,e2)=msh%rhs(:,1,e2)-flx

        ! Left flux boundary
        ql = msh%q(1,msh%nshp,e1)
        qr = msh%q(1,1,i)
        call flux(ql,qr,flx)
        if(i.eq.50) write(*,*) '  flux 1:', ql,qr,flx*msh%detJ(i)
        if (msh%iblank(msh%e2n(2,e1)) > 0) msh%rhs(:,1,i)=msh%rhs(:,1,i)+flx*msh%detJ(i)
        
        ! Right flux boundary
        ql = msh%q(1,msh%nshp,i)
        qr = msh%q(1,1,e2)
        call flux(ql,qr,flx)
        if(i.eq.50) write(*,*) '  flux 2:', ql,qr,flx*msh%detJ(i)
        if (msh%iblank(msh%e2n(1,e2)) > 0) msh%rhs(:,msh%nshp,i)=msh%rhs(:,msh%nshp,i)-flx*msh%detJ(i)
        if(i.eq.50) write(*,*) '  rhs:', msh%rhs(1,1,i),msh%rhs(1,msh%nshp,i)
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
     
     
