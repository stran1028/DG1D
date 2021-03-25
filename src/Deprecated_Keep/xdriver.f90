program conservative_overset
  use code_types
  use pde
  use overset
  implicit none
  !
  integer, parameter :: nmesh=1 
  integer :: i,ntime,n
  real*8 :: dt,mom1(2,nmesh),mom0(2,nmesh)
  real*8, allocatable :: elemInfo1(:),elemInfo2(:)
  integer :: nincomp1,nincomp2
  !
  type(mesh), allocatable :: msh(:)
  allocate(msh(nmesh))
  !
  ntime= 100
  dt=0.25d0 !0.005d0
  !
  call set_type('linear_advection',1d0)
  !call set_type('burgers')
  !
  call init_mesh(msh(1),[-100d0,100d0],1d0,1) ![-1d0,1d0],0.02d0,1)
! call init_mesh(msh(2),[-0.268d0,0.513d0],0.012d0,0)
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
   do n=1,nmesh
      call computeRHS(msh(n))
   enddo
   ! Adjust RHS based on overlaps
   if (nmesh > 1) then
      call findIncompleteElements(msh(1),elemInfo1,nincomp1)
      call findIncompleteElements(msh(2),elemInfo2,nincomp2)
      write(*,*) ' '
      write(*,*) 'ST DEBUG fix incomp 1'
      call fixfluxIncompleteElements(msh(1),msh(2),elemInfo2,nincomp2)
      write(*,*) ' '
      write(*,*) 'ST DEBUG fix incomp 2'
      call fixfluxIncompleteElements(msh(2),msh(1),elemInfo1,nincomp1)
!      call setRHS(msh(1),elemInfo1,nincomp1)
!      call setRHS(msh(2),elemInfo2,nincomp2)
   end if
   do n=1,nmesh
      if (i==1) then
        call computeMoments(msh(n),mom0(:,n))
        !write(6,*) mom0(1,:)
        !write(6,*) mom0(2,:)
        !mom0(1,:)=[1.9165265950903679E-002,  0.37717732514274477 ]
        !mom0(2,:)=[1.2355076923076924,       0.76449230769230769 ]
      endif
      ! Advance time
      call timeIntegrate(msh(n),dt)
      call computeMoments(msh(n),mom1(:,n))
   enddo
   !write(20,*) (mom1(1,1)+mom1(1,2)-mom0(1,1)-mom0(1,2))
   !write(6,*) mom1(2,1)+mom1(2,2)
  end do
  !
  ! write final output
  do n=1,nmesh
     call output(nmesh+n,msh(n))
  enddo
  !
end program conservative_overset
