subroutine init_mesh(msh,xlim,dx,iperiodic)
  !
  use code_types
  use bases
  implicit none
  !
  type(mesh), intent(inout) :: msh
  real*8, intent(in) :: xlim(2)
  real*8, intent(in) :: dx
  integer, intent(in) :: iperiodic
  !
  real*8 :: dxmod
  integer:: i,j,k,ii
  !
  msh%nelem=nint((xlim(2)-xlim(1))/dx)
  dxmod=(xlim(2)-xlim(1))/msh%nelem
  !
  msh%nfields=1
  msh%porder=1
  msh%iperiodic=iperiodic
  msh%nshp = msh%porder+1 
  !
  msh%ngauss=3
  allocate(msh%xgauss(msh%ngauss))
  ! fixed 2nd order gauss pts for now
  !msh%xgauss = [0d0,-sqrt(3d0)/2d0,sqrt(3d0)/2d0] ! Jay's original
  msh%xgauss = 0.5*[0d0,-sqrt(3d0/5d0),sqrt(3d0/5d0)] ! In coordinates of parent element
  msh%wgauss = 0.5*[8d0/9d0,5d0/9d0,5d0/9d0]
!  msh%xgauss = 0.5d0*[0d0,-0.90618d0,-0.538469d0,.538469d0,0.90618d0]
!  msh%wgauss = 0.5d0*[0.568889d0, 0.236927d0,0.478629d0,0.478629d0,0.236927d0]

  msh%nnodes=msh%nelem*2
  allocate(msh%x(msh%nnodes))
  allocate(msh%e2n(2,msh%nelem)) ! element to node number map? 
  allocate(msh%xe(2,msh%nelem))
  allocate(msh%dx(msh%nelem))
  do i=1,msh%nelem
     msh%x(2*i-1) = xlim(1)+(i-1)*dxmod ! left node cood
     msh%x(2*i) = xlim(1)+(i)*dxmod ! right node cood
     msh%e2n(1,i)=2*i-1                 ! left node id
     msh%e2n(2,i)=2*i                   ! right node id
     msh%xe(1,i)=msh%x(2*i-1)               ! left node coord
     msh%xe(2,i)=msh%x(2*i)             ! right node coord
     msh%dx(i) = msh%xe(2,i)-msh%xe(1,i)
  enddo
  !
  allocate(msh%shp(msh%ngauss,msh%nshp))
  allocate(msh%dshp(msh%ngauss,msh%nshp))  
  do i = 1,msh%ngauss
    ! 1D Lagrange elements for now
    ! Parent element goes from -1/2 to 1/2
    call shapefunction(msh%nshp,msh%xgauss(i),[-0.5d0,0.5d0],[1d0,1d0],msh%shp(i,:),msh%dshp(i,:))
  enddo
  !
  msh%nfaces=msh%nelem
  allocate(msh%face(2,msh%nelem))
  !
  do i=1,msh%nelem
     msh%face(1,i)=max(i-1,1) ! gives you left element id
     msh%face(2,i)=min(i+1,msh%nelem) ! gives you right element id
  enddo
  !
  if (iperiodic==1) then
     msh%face(1,1)=msh%nelem
     msh%face(2,msh%nelem)=1
  endif
  !
  allocate(msh%q(msh%nfields,msh%nshp,msh%nelem))
  allocate(msh%rhs(msh%nfields,msh%nshp,msh%nelem))
  allocate(msh%rhsv(msh%nfields,msh%nshp,msh%nelem))
  allocate(msh%rhsf(msh%nfields,msh%nshp,msh%nelem))
  allocate(msh%mass(msh%nfields,msh%nshp,msh%nshp,msh%nelem))
  allocate(msh%iblank(msh%nnodes))
  allocate(msh%nres(msh%nnodes))
  !
  ! Storing NaNb for M = int(Na Nb |J| dx)
  msh%mass = 0.0
  do ii = 1,msh%nelem
  do i = 1,msh%nshp
  do j = 1,msh%nshp
    do k = 1,msh%ngauss
      msh%mass(1,i,j,ii) = msh%mass(1,i,j,ii) + msh%shp(k,i)*msh%shp(k,j)*msh%wgauss(k)*msh%dx(ii)
    enddo
  enddo
  enddo
  enddo
  !
  msh%iblank=1
  msh%nres=dxmod
  !
  if (iperiodic==0) then
     msh%nres(1)=1E15
     msh%nres(msh%nnodes)=1E15
  endif
  !
end subroutine init_mesh
