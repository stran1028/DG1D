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
  real*8 :: err0(nmesh),err1(nmesh)
  real*8, allocatable :: elemInfo1(:),elemInfo2(:)
  integer :: nincomp1,nincomp2,nrk
  real*8 :: rk(4),dx(nmesh),ainf,cfl
  integer :: h
  !
  type(mesh), allocatable :: msh(:)
  allocate(msh(nmesh))
  !
  ! Inputs
  cfl = 0.01d0
  ainf = 1d0
  consoverset = 1

  if(consoverset.eq.1) then 
    write(*,*) 'CONSERVATIVE OVERSET:'
  else
    write(*,*) 'BASELINE (ABUTTING METHOD):'
  endif

  ! Set up the problem and bases types
  call set_type('linear_advection',ainf)
  !call set_type('burgers')
  !call setshp('lagrange')
  call setshp('legendre')
  !
  do order = 0,7
  do h = 1,4
    if(h.eq.1) then  
      dx = [0.2d0,0.1d0]
    elseif(h.eq.2) then 
      dx = [0.08d0,0.04d0]
    elseif(h.eq.3) then 
      dx = [0.01d0,0.005d0]
    elseif(h.eq.4) then 
      dx = [0.02d0,0.01d0]
    endif

    ! Compute parameters
    dt=cfl*minval(dx)/ainf
    ntime = nint(2d0/(ainf*dt)) ! assuming lenght of domain is 2
    write(*,*) '----------------------------------'
    write(*,*) '  INPUTS: '
    write(*,*) '    h, dx = ',h,dx
    write(*,*) '    p = ',order
    write(*,*) '    cfl = ',cfl
    write(*,*) '    dt = ',dt
    write(*,*) '    ntime = ',ntime
    write(*,*) '    test = ',cfl,minval(dx),ainf
    !
    ! Initialize the mesh(es)
    call init_mesh(msh(1),[-1d0,1d0],dx(1),1,order)
    call init_mesh(msh(2),[-0.268d0,0.732d0],dx(2),0,order)
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
        call fixMassIncompleteElements(msh(1),msh(2),elemInfo2,nincomp2)
        call fixMassIncompleteElements(msh(2),msh(1),elemInfo1,nincomp1)
      endif
    end if
    !
    do n=1,nmesh
     !call output(n,msh(n))
     call output(100*order+10*h+n,msh(n))
     call computeMoments(msh(n),mom0(:,n),err0(n))
    enddo
    !
    ! Iterate in time
    rk = [1d0/4d0, 8d0/15d0,5d0/12d0, 3d0/4d0];
    do i=1,ntime
   
!     write(*,*) '--------------------------'
!     write(*,*) 'TIMESTEP ',i
!     write(*,*) '--------------------------'
     ! RK step 1
     call timestep(nmesh,dt,msh,consoverset,elemInfo1,elemInfo2,nincomp1,nincomp2)
     do j = 1,nmesh
       ! Euler 1st order
!       msh(j)%q=msh(j)%q+dt*msh(j)%dq
!       msh(j)%sol=msh(j)%q

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
!       call output(nmesh+n,msh(n))
       call output(100*order+10*h+nmesh+n,msh(n))
       call computeMoments(msh(n),mom1(:,n),err1(n))
    enddo
    write(*,*) '  Initial L2 error: ',sum(err0(:))
    write(*,*) '  Conservation error: ',sum(mom0(1,:)),sum(mom1(1,:)),sum(mom1(1,:))-sum(mom0(1,:))
    write(*,*) '  Final L2 error, diff: ',sum(err1(:)),sum(err1(:))-sum(err0(:))
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
end program conservative_overset
