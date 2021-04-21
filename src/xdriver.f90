program conservative_overset
  use code_types
  use pde
  use overset
  use bases
  implicit none
  !
  integer, parameter :: nmesh=2
  integer :: i,ntime,n,j,k,m,order,consoverset
  real*8 :: dt,mom1(2,nmesh),mom0(2,nmesh)
  real*8, allocatable :: elemInfo1(:),elemInfo2(:)
  integer :: nincomp1,nincomp2,nrk
  real*8 :: rk(4),dx(nmesh),ainf,cfl
  !
  type(mesh), allocatable :: msh(:)
  allocate(msh(nmesh))
  !
  ! Inputs
  cfl = 0.1d0
  dx = [0.02d0,0.01d0]
  ainf = 1
  order = 2
  consoverset = 0

  ! Compute parameters
  dt=cfl*minval(dx)/ainf
  ntime = nint(2d0/(ainf*dt)) ! assuming lenght of domain is 2
  write(*,*) 'INPUTS: '
  write(*,*) '  cfl = ',cfl
  write(*,*) '  dt = ',dt
  write(*,*) '  dx = ',dx
  write(*,*) '  ntime = ',ntime
  !
  ! Set up the problem and bases types
  call set_type('linear_advection',ainf)
  !call set_type('burgers')
!  call setshp('lagrange')
  call setshp('legendre')
  !
  ! Initialize the mesh(es)
  call init_mesh(msh(1),[-1d0,1d0],dx(1),1,order)
  call init_mesh(msh(2),[-0.268d0,0.518d0],dx(2),0,order)
  !
  do n=1,nmesh
   call initvar(msh(n))
  enddo
  !
  ! Blank out coarser overlapping cells 
  if(nmesh>1) then 
    call connect(msh(1),msh(2))
    call connect(msh(2),msh(1))
    call fixOverlap(msh(2),msh(1))
    call findIncompleteElements(msh(1),elemInfo1,nincomp1)
    call findIncompleteElements(msh(2),elemInfo2,nincomp2)
    if(consoverset.eq.1) then 
    write(*,*) ' '
      call fixMassIncompleteElements(msh(1),msh(2),elemInfo2,nincomp2)
      call fixMassIncompleteElements(msh(2),msh(1),elemInfo1,nincomp1)
    endif
  end if
  !
  ! write IC
  write(*,*) '----------------------------------'
  write(*,*) "INITIAL CONDITION:"
  do n=1,nmesh
   call output(n,msh(n))
   call computeMoments(msh(n),mom0(:,n),ainf,dt,ntime)
  enddo
  write(*,*) '----------------------------------'
  !
  ! Iterate in time
  rk = [1d0/4d0, 8d0/15d0,5d0/12d0, 3d0/4d0];
  do i=1,ntime
   ! Compute RHS for all elements
!   write(*,*) ' '
!   write(*,*) '==================='
!   write(*,*) 'Step ',i
!   write(*,*) '==================='
   
   ! RK step 1
   call timestep(nmesh,dt,msh,consoverset,elemInfo1,elemInfo2,nincomp1,nincomp2)
   do j = 1,nmesh
     ! Euler 1st order
!     msh(j)%q=msh(j)%q+dt*msh(j)%dq
!     msh(j)%sol=msh(j)%q

     msh(j)%q=msh(j)%sol+rk(2)*dt*msh(j)%dq
     msh(j)%sol=msh(j)%sol+rk(1)*dt*msh(j)%dq
   enddo

   ! RK step 2
   call timestep(nmesh,dt,msh,consoverset,elemInfo1,elemInfo2,nincomp1,nincomp2)
   do j = 1,nmesh
     msh(j)%q=msh(j)%sol+rk(3)*dt*msh(j)%dq
   enddo

  ! RK step 3
   call timestep(nmesh,dt,msh,consoverset,elemInfo1,elemInfo2,nincomp1,nincomp2)
   do j = 1,nmesh
     msh(j)%sol=msh(j)%sol+rk(4)*dt*msh(j)%dq
     msh(j)%q = msh(j)%sol
   enddo

   ! check solution
  ! do j = 1,nmesh
  ! do m = 1,msh(j)%nelem
  ! do k = 1,msh(j)%nshp
  !   if(isnan(msh(j)%q(1,k,m))) then 
  !     write(*,*) 'Found Nan at ',i,j,k,m
  !     call exit(1)
  !   endif
  ! enddo
  ! enddo
  ! enddo

  end do ! timesteps
  !
  ! write final output
  write(*,*) '----------------------------------'
  write(*,*) 'FINAL OUTPUT: '
  do n=1,nmesh
     call output(nmesh+n,msh(n))
     call computeMoments(msh(n),mom1(:,n))
  enddo
  !
end program conservative_overset
