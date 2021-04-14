program conservative_overset
  use code_types
  use pde
  use overset
  use bases
  implicit none
  !
  integer, parameter :: nmesh=2
  integer :: i,ntime,n,j,order
  real*8 :: dt,mom1(2,nmesh),mom0(2,nmesh)
  real*8, allocatable :: elemInfo1(:),elemInfo2(:)
  integer :: nincomp1,nincomp2,nrk
  real*8 :: rk(4)
  !
  type(mesh), allocatable :: msh(:)
  allocate(msh(nmesh))
  !
  ntime= 4000 ! 4000
  dt=.1d0 !0.010471975511966d0 !0.05d0/3d0
  rk = [1d0/4d0, 8d0/15d0,5d0/12d0, 3d0/4d0];
  !
  ! Set up the problem and bases types
  call set_type('linear_advection',1d0)
  !call set_type('burgers')
!  call setshp('lagrange')
  call setshp('legendre')
  !
  ! Initialize the mesh(es)
!call init_mesh(msh(1),[0d0,6.28318530717959d0],0.314159265358979d0,1)
!call init_mesh(msh(1),[-1d0,1d0],0.05d0,1)

  order = 1
  call init_mesh(msh(1),[-100d0,100d0],1d0,1,order)
  call init_mesh(msh(2),[-10.1d0,20.6d0],0.5d0,0,order)
! call init_mesh(msh(2),[-15.5d0,50.5d0],2d0,0)
  !
  do n=1,nmesh
   call initvar(msh(n))
  enddo
  !
  ! Blank out coarser overlapping cells 
  if(nmesh>1) then 
    write(*,*) 'Connect 1 and 2'
    call connect(msh(1),msh(2))
    write(*,*) 'Connect 2 and 1'
    call connect(msh(2),msh(1))
    call fixOverlap(msh(2),msh(1))
    call findIncompleteElements(msh(1),elemInfo1,nincomp1)
    call findIncompleteElements(msh(2),elemInfo2,nincomp2)
    call fixMassIncompleteElements(msh(1),msh(2),elemInfo2,nincomp2)
    call fixMassIncompleteElements(msh(2),msh(1),elemInfo1,nincomp1)
  end if
  !
  ! write IC
  do n=1,nmesh
   call output(n,msh(n))
  enddo
  !
  ! Iterate in time
  do i=1,ntime
   ! Compute RHS for all elements
   write(*,*) ' '
   write(*,*) '==================='
   write(*,*) 'Step ',i
   write(*,*) '==================='
   
   ! RK step 1
   call timestep(nmesh,dt,msh)
   do j = 1,nmesh
     msh(j)%q=msh(j)%sol+rk(2)*dt*msh(j)%dq
     msh(j)%sol=msh(j)%sol+rk(1)*dt*msh(j)%dq
   enddo

   ! RK step 2
   call timestep(nmesh,dt,msh)
   do j = 1,nmesh
     msh(j)%q=msh(j)%sol+rk(3)*dt*msh(j)%dq
   enddo

   ! RK step 3
   call timestep(nmesh,dt,msh)
   do j = 1,nmesh
     msh(j)%sol=msh(j)%sol+rk(4)*dt*msh(j)%dq
     msh(j)%q = msh(j)%sol
   enddo

  end do ! timesteps
  !
  ! write final output
  do n=1,nmesh
     call output(nmesh+n,msh(n))
  enddo
  !
end program conservative_overset
