program conservative_overset
  use code_types
  use pde
  use overset
  use bases
  implicit none
  !
  integer, parameter :: nmesh=1
  integer :: i,ntime,n,j
  real*8 :: dt,mom1(2,nmesh),mom0(2,nmesh)
  real*8, allocatable :: elemInfo1(:),elemInfo2(:)
  integer :: nincomp1,nincomp2,nrk
  real*8 :: rk(4)
  !
  type(mesh), allocatable :: msh(:)
  allocate(msh(nmesh))
  !
  ntime= 20000
  dt=0.010471975511966d0 !0.05d0/3d0
  nrk = 3
  rk = [1d0/4d0, 8d0/15d0,5d0/12d0, 3d0/4d0];
  write(*,*) 'rk = ',rk(2)
  !
  ! Set up the problem and bases types
  call set_type('linear_advection',1d0)
  !call set_type('burgers')
  call setshp('lagrange')
!  call setshp('legendre')
  !
  ! Initialize the mesh(es)
call init_mesh(msh(1),[0d0,6.28318530717959d0],0.314159265358979d0,1)
!call init_mesh(msh(1),[-1d0,1d0],0.05d0,1)

!  call init_mesh(msh(1),[-100d0,100d0],2d0,1)
!  call init_mesh(msh(2),[-25.5d0,50.5d0],1d0,0)
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
  endif
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
   write(*,*) 'q0 = ', msh(1)%q(1,:,10)
   write(*,*) ' '
   call timestep(nmesh,dt,msh)
   do j = 1,nmesh
     msh(j)%q=msh(j)%sol+rk(2)*dt*msh(j)%dq
     msh(j)%sol=msh(j)%sol+rk(1)*dt*msh(j)%dq
   enddo
   write(*,*) 'dq1 = ', msh(1)%dq(1,:,10)
   write(*,*) 'q1 = ', msh(1)%q(1,:,10)
   write(*,*) 'sol1 = ', msh(1)%sol(1,:,10)
   write(*,*) ' '

   ! RK step 2
   call timestep(nmesh,dt,msh)
   do j = 1,nmesh
     msh(j)%q=msh(j)%sol+rk(3)*dt*msh(j)%dq
!     msh(j)%sol = msh(j)%q
   enddo
   write(*,*) 'dq2 = ', msh(1)%dq(1,:,10)
   write(*,*) 'q2 = ', msh(1)%q(1,:,10)
   write(*,*) 'sol2 = ', msh(1)%sol(1,:,10)
   write(*,*) ' '

   ! RK step 3
   call timestep(nmesh,dt,msh)
   do j = 1,nmesh
     msh(j)%sol=msh(j)%sol+rk(4)*dt*msh(j)%dq
   enddo
   write(*,*) 'dq3 = ', msh(1)%dq(1,:,10)
   write(*,*) 'q3 = ', msh(1)%q(1,:,10)
   write(*,*) 'sol3 = ', msh(1)%sol(1,:,10)
   write(*,*) ' '
  end do ! timesteps
  !
  ! write final output
  do n=1,nmesh
     call output(nmesh+n,msh(n))
  enddo
  !
end program conservative_overset
