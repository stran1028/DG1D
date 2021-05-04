program conservative_overset
  use code_types
  use pde
  use overset
  use bases
  implicit none
  !
  integer, parameter :: nmesh=2
  integer :: i,s,ntime,n,j,k,m,order,consoverset
  real*8 :: dt,mom1(2,nmesh),mom0(2,nmesh)
  real*8 :: err(nmesh)
  real*8, allocatable :: elemInfo1(:,:),elemInfo2(:,:)
  integer :: nincomp1,nincomp2,nrk
  integer :: conswitch,noverlap
  real*8 :: rk(4),dx(nmesh),ainf,cfl,foverlap
  integer :: h
  !
  type(mesh), allocatable :: msh(:)
  allocate(msh(nmesh))
  !
  ! Inputs
  cfl = 0.01d0
  ainf = 1d0
  foverlap = 0.0 ! 0 = coarse mesh clips all, 0.5 = both clip half, 1 = fine mesh clip all
  !
  if((foverlap.gt.1d0).or.(foverlap.lt.0d0)) then 
    write(*,*) 'foverlap wrong. try again.'
    call exit(1)
  endif
  !

  ! Set up the problem and bases types
  call set_type('linear_advection',ainf)
  !call set_type('burgers')
  !
  do conswitch = 1,1
  do s = 1,2
  do noverlap = 3,3
  do order = 1,5
  do h = 0,5

    if ((conswitch.eq.0).and.(noverlap.gt.1)) cycle 

    write(*,*) '----------------------------------'
    write(*,*) 'INPUTS: '

    if(conswitch.eq.0) then 
      consoverset = 0
      write(*,*) '    BASELINE (ABUTTING METHOD):'
    else
      consoverset = 1
      write(*,*) '    CONSERVATIVE OVERSET:'
    endif 

    if(noverlap.eq.1) then 
       foverlap = 0.0d0         ! Coarse mesh clips all of overlap
    elseif(noverlap.eq.2) then 
       foverlap = 0.5d0         ! Both meshes clip 0.5 overlap each
    elseif(noverlap.eq.3) then 
       foverlap = 1.0d0         ! Fine mesh clips all of overlap
    endif
    write(*,*) '    foverlap = ',foverlap
    
    if(s.eq.1) then 
      call setshp('lagrange')
      write(*,*) '    LAGRANGE SHAPE FUNCTIONS'
    else
      call setshp('legendre')
      write(*,*) '    LEGENDRE SHAPE FUNCTIONS'
    endif

    if(h.eq.0) then  
      dx = [0.5d0,0.25d0]
    elseif(h.eq.1) then  
      dx = [0.25d0,0.125d0]
    elseif(h.eq.2) then 
      dx = [0.125d0,0.0625d0]
    elseif(h.eq.3) then 
      dx = [0.0625d0,0.03125d0]
    elseif(h.eq.4) then 
      dx = [0.03125d0,0.015625d0]
    elseif(h.eq.5) then 
      dx = [0.015625d0,0.0078125d0]
    endif

    ! Compute parameters
    dt=cfl*minval(dx)/ainf
    ntime =   nint(2d0/(ainf*dt)) ! assuming lenght of domain is 2
    write(*,*) '    h, dx = ',h,dx
    write(*,*) '    p = ',order
    write(*,*) '    cfl = ',cfl
    write(*,*) '    dt = ',dt
    write(*,*) '    ntime = ',ntime
    !
    ! Initialize the mesh(es)
    call init_mesh(msh(1),[-1d0,1d0],dx(1),1,order)
    call init_mesh(msh(2),[-0.268d0,0.732d0],dx(2),0,order)
    !
    do n=1,nmesh
     call initvar(msh(n))
     msh(n)%sol=msh(n)%q
    enddo
    ! Store initial conditions (exact solution)
    msh(1)%q0 = msh(1)%q 
    msh(2)%q0 = msh(2)%q 

    !
    ! Blank out coarser overlapping cells 
    if(nmesh>1) then 
      call connect(msh(1),msh(2))
      call connect(msh(2),msh(1))
      call fixOverlap(msh(2),msh(1))
      call findIncompleteElements(msh(1),elemInfo1,nincomp1)
      call findIncompleteElements(msh(2),elemInfo2,nincomp2)
      call fixMassIncompleteElements(msh(1),msh(2),elemInfo2,nincomp2,&
              consoverset,foverlap)
      call fixMassIncompleteElements(msh(2),msh(1),elemInfo1,nincomp1,&
              consoverset,foverlap)
    end if
    !

    do n=1,nmesh
!     call output(100*order+10*h+n,msh(n))
    enddo
    call computeMoments(msh(1),mom0(:,1),err(1),nincomp1,elemInfo1)
    call computeMoments(msh(2),mom0(:,2),err(2),nincomp2,elemInfo2)
    !
    ! Iterate in time
    rk = [1d0/4d0, 8d0/15d0,5d0/12d0, 3d0/4d0];
    do i=1,ntime
   
!     write(*,*) '--------------------------'
!     write(*,*) 'TIMESTEP ',i
!     write(*,*) '--------------------------'
     ! RK step 1
     call timestep(nmesh,dt,msh,consoverset,elemInfo1,elemInfo2,nincomp1,nincomp2,foverlap)

     do j = 1,nmesh
       ! Euler 1st order
!       msh(j)%q=msh(j)%q+dt*msh(j)%dq
!       msh(j)%sol=msh(j)%q
      msh(j)%q=msh(j)%sol+rk(2)*dt*msh(j)%dq
      msh(j)%sol=msh(j)%sol+rk(1)*dt*msh(j)%dq
     enddo

     ! RK step 2
     call timestep(nmesh,dt,msh,consoverset,elemInfo1,elemInfo2,nincomp1,nincomp2,foverlap)
     do j = 1,nmesh
      msh(j)%q=msh(j)%sol+rk(3)*dt*msh(j)%dq
     enddo

    ! RK step 3
     call timestep(nmesh,dt,msh,consoverset,elemInfo1,elemInfo2,nincomp1,nincomp2,foverlap)
     do j = 1,nmesh
      msh(j)%sol=msh(j)%sol+rk(4)*dt*msh(j)%dq
      msh(j)%q = msh(j)%sol
     enddo

  !   check solution
  !   do j = 1,nmesh
  !   do m = 1,msh(j)%nelem
  !   do k = 1,msh(j)%nshp
  !      if(isnan(msh(j)%q(1,k,m))) then 
  !       write(*,*) 'Found Nan at ',i,j,k,m
  !       call exit(1)
  !     endif
  !   enddo
  !   enddo
  !   enddo

    end do ! timesteps
    !
    ! write final output
    write(*,*) ' '
    write(*,*) '  FINAL OUTPUT: '
    do n=1,nmesh
!       call output(100*order+10*h+nmesh+n,msh(n))
    enddo
    call computeMoments(msh(1),mom1(:,1),err(1),nincomp1,elemInfo1)
    call computeMoments(msh(2),mom1(:,2),err(2),nincomp2,elemInfo2)

    write(*,*) '  Conservation error: ',sum(mom0(1,:)),sum(mom1(1,:)),sum(mom1(1,:))-sum(mom0(1,:))
    write(*,*) '  Final L2 error: ',sqrt(sum(err(:)))
    write(*,*) '----------------------------------'
    !
    ! Clear memory
    do n = 1,nmesh
      call clearMem(msh(n))
    enddo
    if (allocated(elemInfo1)) deallocate(elemInfo1)
    if (allocated(elemInfo2)) deallocate(elemInfo2)

  enddo ! mesh sweep
  enddo ! p sweep
  enddo
  enddo
  enddo
end program conservative_overset
