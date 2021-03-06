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
  real*8 :: rk(4),dx(nmesh),ainf,cfl,foverlap,sweep(5,2)
  real*8 :: test1(6),test2(6)
  real*8 :: time(2),m2start
  integer :: h
  !
  type(mesh), allocatable :: msh(:)
  allocate(msh(nmesh))
  !
  ! Inputs
  cfl = 0.01d0
  ainf = 1d0
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
  do conswitch = 0,1    ! cons overset loop 
  do s = 2,2            ! shape function loop
  do noverlap = 1,3     ! foverlap loop
  do order = 1,5        ! p-order loop
    sweep = 0d0

    if ((conswitch.eq.0).and.(noverlap.gt.1)) cycle 

    write(*,*) '----------------------------------'
    write(*,*) 'INPUTS: '

    if(conswitch.eq.0) then 
      consoverset = 0
      write(*,*) '  BASELINE (ABUTTING METHOD):'
    else
      consoverset = 1
      write(*,*) '  CONSERVATIVE OVERSET:'
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
      write(*,*) '  LAGRANGE SHAPE FUNCTIONS'
    else
      call setshp('legendre')
      write(*,*) '  LEGENDRE SHAPE FUNCTIONS'
    endif

    ! Do a mesh sweep
    do h = 1,5
      ! start timer
      call cpu_time(time(1))
      if(h.eq.1) then  
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

      ! DEBUG, do an overlap sweep with a constant mesh size
!      dx = [0.25d0,0.125d0]
!      if(h.eq.1) then  
!        m2start = -0.5 - dx(2)*.25
!      elseif(h.eq.2) then  
!        m2start = -0.5 - dx(2)*.50
!      elseif(h.eq.3) then 
!        m2start = -0.5 - dx(2)*.75
!      elseif(h.eq.4) then 
!        m2start = -0.5 - dx(2)*.90
!      elseif(h.eq.5) then 
!        m2start = -0.5 - dx(2)*.95
!      elseif(h.eq.6) then 
!        m2start = -0.5 - dx(2)*.99
!      endif

      m2start = -0.5 - dx(2)*.95 !! stress test w/ 95% cut

      ! Compute parameters
      dt=cfl*minval(dx)/ainf
      ntime =  nint(2d0/(ainf*dt)) ! assuming lenght of domain is 2
      write(*,*) ' '
      write(*,*) '    h, dx = ',h,dx
      write(*,*) '    m2start = ',m2start
      write(*,*) '    p = ',order
      write(*,*) '    cfl = ',cfl
      write(*,*) '    dt = ',dt
      write(*,*) '    ntime = ',ntime
      !
      ! Initialize the mesh(es)
      call init_mesh(msh(1),[-1d0,1d0],dx(1),1,order)
      call init_mesh(msh(2),[m2start,m2start+1d0],dx(2),0,order)
!      call init_mesh(msh(2),[-0.268d0,0.732d0],dx(2),0,order)
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
       call output(n,msh(n))
      enddo
      call computeMoments(msh(1),mom0(:,1),err(1),nincomp1,elemInfo1)
      call computeMoments(msh(2),mom0(:,2),err(2),nincomp2,elemInfo2)
      !
      ! Iterate in time
      rk = [1d0/4d0, 8d0/15d0,5d0/12d0, 3d0/4d0];
      do i=1,ntime
       !write(*,*) '--------------------------'
       !write(*,*) 'TIMESTEP ',i
       !write(*,*) '--------------------------'
       ! RK step 1
        call timestep(nmesh,dt,msh,consoverset,elemInfo1,elemInfo2,nincomp1,nincomp2,foverlap)

        do j = 1,nmesh
          ! Euler 1st order
!         msh(j)%q=msh(j)%q+dt*msh(j)%dq
!         msh(j)%sol=msh(j)%q

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
       enddo ! timesteps
       !
       ! write final output
       do n=1,nmesh
        call output(nmesh+n,msh(n))
       enddo
       call computeMoments(msh(1),mom1(:,1),err(1),nincomp1,elemInfo1)
       call computeMoments(msh(2),mom1(:,2),err(2),nincomp2,elemInfo2)
       sweep(h,1) = sqrt(sum(err(:)))
       sweep(h,2) = sum(mom1(1,:))-sum(mom0(1,:))
       write(*,*) '    Min Rem Frac M1: ',minval(msh(1)%dxcut)/dx(1)
       write(*,*) '    Min Rem Frac M2: ',minval(msh(2)%dxcut)/dx(2)
       !
       ! Clear memory
       do n = 1,nmesh
         call clearMem(msh(n))
       enddo
       if (allocated(elemInfo1)) deallocate(elemInfo1)
       if (allocated(elemInfo2)) deallocate(elemInfo2)

       ! call timer
       call cpu_time(time(2))
       write(*,*) ' '
       write(*,*) '    Run time: ',time(2)-time(1)
       write(*,*) '    Run time/step: ',(time(2)-time(1))/ntime
     enddo ! mesh sweep

     write(*,*) ' '
     write(*,*) 'FINAL SWEEP OUTPUT: '
     write(*,*) '  H SWEEP L2 ERROR = ',sweep(:,1)
     write(*,*) '  H SWEEP CONS ERROR = ',sweep(:,2)
     write(*,*) '----------------------------------'
  enddo ! p sweep
  enddo
  enddo
  enddo
end program conservative_overset
