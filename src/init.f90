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
  msh%ngauss= msh%nshp+1
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
  elseif(msh%ngauss.eq.5) then 
    msh%xgauss = 0.5*[0.0000000,-0.5384693101056831,0.5384693101056831,-0.9061798459386640,0.9061798459386640]
    msh%wgauss = 0.5*[0.5688888888888889,0.4786286704993665,0.4786286704993665,0.2369268850561891,0.2369268850561891]
  elseif(msh%ngauss.eq.6) then 
    msh%xgauss = 0.5*[-0.6612093864662645,0.6612093864662645,-0.2386191860831969,0.2386191860831969,&
                      -0.9324695142031521,0.9324695142031521]
    msh%wgauss = 0.5*[0.3607615730481386,0.3607615730481386,0.4679139345726910,0.4679139345726910,&
                      0.1713244923791704,0.1713244923791704]
  elseif(msh%ngauss.eq.7) then 
    msh%xgauss = 0.5*[0.0000000,-0.4058451513773972,0.4058451513773972,-0.7415311855993945,0.7415311855993945,&
                     -0.9491079123427585,0.9491079123427585]
    msh%wgauss = 0.5*[0.4179591836734694,0.3818300505051189,0.3818300505051189,0.2797053914892766,&
                      0.2797053914892766,0.1294849661688697,0.1294849661688697]
  elseif(msh%ngauss.eq.8) then 
    msh%xgauss = 0.5*[-0.1834346424956498,0.1834346424956498,-0.5255324099163290,0.5255324099163290,&
                      -0.7966664774136267,0.7966664774136267,-0.9602898564975363,0.9602898564975363]
    msh%wgauss = 0.5*[0.3626837833783620,0.3626837833783620,0.3137066458778873,0.3137066458778873,&
                      0.2223810344533745,0.2223810344533745,0.1012285362903763,0.1012285362903763]
  elseif(msh%ngauss.eq.9) then 
    msh%xgauss = 0.5*[0.000000,-0.8360311073266358,0.8360311073266358,-0.9681602395076261,0.9681602395076261,&
                     -0.3242534234038089,0.3242534234038089,-0.6133714327005904,0.6133714327005904]
    msh%wgauss = 0.5*[0.3302393550012598,0.1806481606948574,0.1806481606948574,0.0812743883615744,&
                      0.0812743883615744,0.3123470770400029,0.3123470770400029,0.2606106964029354,0.2606106964029354]
  elseif(msh%ngauss.eq.10) then 
    msh%xgauss = 0.5*[-0.1488743389816312,0.1488743389816312,-0.4333953941292472,0.4333953941292472,&
                      -0.6794095682990244,0.6794095682990244,-0.8650633666889845,0.8650633666889845,&
                      -0.9739065285171717,0.9739065285171717]
    msh%wgauss = 0.5*[0.2955242247147529,0.2955242247147529,0.2692667193099963,0.2692667193099963,&
                      0.2190863625159820,0.2190863625159820,0.1494513491505806,0.1494513491505806,&
                      0.0666713443086881,0.0666713443086881]
  endif
  write(*,*) '    order,nelem,nshp,ngauss = ',order,msh%nelem,msh%nshp,msh%ngauss

  msh%nnodes=msh%nelem*msh%nshp
  msh%nvtx=msh%nelem*2
  allocate(msh%x(msh%nnodes))
  allocate(msh%e2n(msh%nshp,msh%nelem)) ! element to node number map? 
  allocate(msh%xe(2,msh%nelem))
  allocate(msh%dx(msh%nelem))

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
       do j = 1,msh%nshp
         tmp = dxmod/(msh%nshp-1d0)
         msh%x(index1+j) = msh%xe(1,i)+tmp*(j-1d0)
       enddo 
     endif 
  enddo
  !
  allocate(msh%shp(msh%ngauss,msh%nshp))
  allocate(msh%dshp(msh%ngauss,msh%nshp))  
  do i = 1,msh%ngauss
    ! 1D Parent element goes from -1/2 to 1/2
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
