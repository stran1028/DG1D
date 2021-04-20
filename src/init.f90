subroutine init_mesh(msh,xlim,dx,iperiodic,order)
  !
  use code_types
  use bases
  implicit none
  !
  type(mesh), intent(inout) :: msh
  real*8, intent(in) :: xlim(2)
  real*8, intent(in) :: dx
  integer, intent(in) :: iperiodic,order
  !
  real*8 :: dxmod,tmp
  integer:: i,j,k,ii,index1
  !
  msh%nelem=nint((xlim(2)-xlim(1))/dx)
  dxmod=(xlim(2)-xlim(1))/msh%nelem
  !
  msh%nfields=1
  msh%porder=order
  msh%iperiodic=iperiodic
  msh%nshp = order+1 
  !
  msh%ngauss=msh%nshp+1
  allocate(msh%xgauss(msh%ngauss))
  allocate(msh%wgauss(msh%ngauss))
  ! Parent element goes from -0.5 to 0.5
  if(msh%ngauss.eq.1) then 
    msh%xgauss = 0d0
    msh%wgauss = 1d0
  elseif(msh%ngauss.eq.2) then 
    msh%xgauss = [-0.5/sqrt(3d0),0.5/sqrt(3d0)]
    msh%wgauss = [0.5,0.5]
  elseif(msh%ngauss.eq.3) then 
    msh%xgauss = 0.5*[0d0,-sqrt(3d0/5d0),sqrt(3d0/5d0)]
    msh%wgauss = 0.5*[8d0/9d0,5d0/9d0,5d0/9d0]
  elseif(msh%ngauss.eq.4) then 
    msh%xgauss = 0.5*[-sqrt(3d0/7d0-2d0/7d0*sqrt(6d0/5d0)),sqrt(3d0/7d0-2d0/7d0*sqrt(6d0/5d0)), &
                      -sqrt(3d0/7d0+2d0/7d0*sqrt(6d0/5d0)),sqrt(3d0/7d0+2d0/7d0*sqrt(6d0/5d0))]
    msh%wgauss = 0.5*[(18d0+sqrt(30d0))/36d0,(18d0+sqrt(30d0))/36d0,(18d0-sqrt(30d0))/36d0,(18d0-sqrt(30d0))/36d0]
  endif
  write(*,*) 'order,nshp,ngauss = ',order,msh%nshp,msh%ngauss
  write(*,*) 'xgauss: ',msh%xgauss
  write(*,*) 'wgauss: ',msh%wgauss

  msh%nnodes=msh%nelem*msh%nshp
  msh%nvtx=msh%nelem*2
  allocate(msh%x(msh%nnodes))
  allocate(msh%e2n(msh%nshp,msh%nelem)) ! element to node number map? 
  allocate(msh%xe(2,msh%nelem))
  allocate(msh%dx(msh%nelem))

! Need to rewrite for higher order elements
  do i=1,msh%nelem
     index1 = (i-1)*msh%nshp
     msh%xe(1,i)=xlim(1)+(i-1)*dxmod        ! left node coord
     msh%xe(2,i)=xlim(1)+i*dxmod        ! right node coord
     msh%dx(i) = msh%xe(2,i)-msh%xe(1,i)
     do j = 1,msh%nshp
       msh%e2n(j,i) = index1 + j ! ids node j on element j
     enddo 

     ! node locations
     if(order.eq.0) then
       msh%x(index1+1) = 0.5*(msh%xe(1,i) + msh%xe(2,i))
     elseif(order.gt.0) then 
       msh%x(index1+1) = msh%xe(1,i)
       do j = 2,msh%nshp
         tmp = dxmod/(msh%nshp-1d0)
         msh%x(index1+j) = msh%xe(1,i)+tmp
       enddo 
     endif 
       
  enddo
  !
  allocate(msh%shp(msh%ngauss,msh%nshp))
  allocate(msh%dshp(msh%ngauss,msh%nshp))  
  do i = 1,msh%ngauss
    ! 1D Lagrange elements for now
    ! Parent element goes from -1/2 to 1/2
    msh%shp(i,:) = 1d0
    call shapefunction(msh%nshp,msh%xgauss(i),[-0.5d0,0.5d0],msh%shp(i,:),msh%shp(i,:),msh%dshp(i,:))
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
  allocate(msh%dq(msh%nfields,msh%nshp,msh%nelem))
  allocate(msh%q(msh%nfields,msh%nshp,msh%nelem))
  allocate(msh%sol(msh%nfields,msh%nshp,msh%nelem))
  allocate(msh%rhs(msh%nfields,msh%nshp,msh%nelem))
  allocate(msh%rhsF(msh%nfields,msh%nshp,msh%nelem))
  allocate(msh%rhsV(msh%nfields,msh%nshp,msh%nelem))
  allocate(msh%mass(msh%nfields,msh%nshp*msh%nshp,msh%nelem))
  allocate(msh%iblank(2,msh%nelem))
  allocate(msh%nres(msh%nvtx))
  !
  ! Storing NaNb for M = int(Na Nb |J| dx)
  msh%mass = 0.0
  do ii = 1,msh%nelem
  do i = 1,msh%nshp
  do j = 1,msh%nshp
    index1 = (i-1)*msh%nshp+j
    do k = 1,msh%ngauss
      msh%mass(1,index1,ii) = msh%mass(1,index1,ii) + msh%shp(k,i)*msh%shp(k,j)*msh%wgauss(k)*msh%dx(ii)
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
     msh%nres(msh%nvtx)=1E15
  endif
  !
end subroutine init_mesh
