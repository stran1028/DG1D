program conservative_overset
  use code_types
  use pde
  use overset
  use bases
  use slopelimiter
  implicit none
  !
  integer, parameter :: nmesh=1
  integer :: i,s,ntime,n,j,k,m,order,consoverset,ilim,ieuler, isupg, ireg
  real*8 :: dt,mom1(2,nmesh),mom0(2,nmesh)
  real*8 :: err(nmesh)
  real*8, allocatable :: elemInfo1(:,:),elemInfo2(:,:)
  integer :: nincomp1,nincomp2,nrk,debug
  integer :: conswitch,noverlap
  real*8 :: rk(4),dx(nmesh),ainf,muinf,cfl,foverlap,sweep(5,3)
  real*8 :: test1(6),test2(6)
  real*8 :: time(2),m2start,q1,qN
  integer :: h,isMerged(nmesh)
  !
  type(mesh), allocatable :: msh(:)
  allocate(msh(nmesh))
  !
  ! Inputs
  cfl = 0.01d0
  ainf = 2d0
  muinf = 0.02d0
  q1 = 2d0
  qN = 0d0
  !
  ! Set up the problem and bases types
!  call set_type('linear_advection',ainf,muinf,q1,qN)
  call set_type('burgers',ainf,muinf,q1,qN)
  ilim = 0      ! flag to control slope limiting
  isupg = 0  ! supg flag
  ireg = 0 ! regularization flag
  ieuler = 0
  do conswitch = 0,0    ! cons overset loop 
  do s = 1,1            ! shape function loop
  do noverlap = 1,1     ! foverlap loop
  do order = 4,4      ! p-order loop
    sweep = 0d0

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
    if ((conswitch.eq.0).and.(noverlap.gt.1)) foverlap = 0.0d0
    write(*,*) '    foverlap = ',foverlap
    
    if(s.eq.1) then 
      call setshp('lagrange')
      write(*,*) '  LAGRANGE SHAPE FUNCTIONS'
    else
      call setshp('legendre')
      write(*,*) '  LEGENDRE SHAPE FUNCTIONS'
    endif

    ! Do a mesh sweep
    do h = 2,2
      ! start timer
      call cpu_time(time(1))
      if(h.eq.1) then  
!        dx = [0.25d0,0.125d0]
      elseif(h.eq.2) then 
!        dx = [0.125d0,0.0625d0]
      elseif(h.eq.3) then 
!        dx = [0.0625d0,0.03125d0]
      elseif(h.eq.4) then 
!        dx = [0.03125d0,0.015625d0]
      elseif(h.eq.5) then 
!        dx = [0.015625d0,0.0078125d0]
      endif
!      m2start = -0.25 - dx(2)*.99 !! stress test w/ 90% cut
     dx = 0.125d0

      ! Compute parameters
      dt=cfl*minval(dx)/ainf
      ntime = 60
!      ntime = nint(1.25d0/dt) ! T = 1.25 seconds for Burgers
!      ntime = 1.*nint(2d0/(ainf*dt)) ! assuming lenght of domain is 2
      write(*,*) ' '
      write(*,*) '    ainf, muinf = ',ainf,muinf
      if (index(pde_descriptor,'burgers') > 0) write(*,*) '    q1, qN = ',q1,qN  
      write(*,*) '    h, dx = ',h,dx
      write(*,*) '    m2start = ',m2start
      write(*,*) '    ilim,isupg,ieuler = ',ilim,isupg,ieuler
      write(*,*) '    p = ',order
      write(*,*) '    cfl = ',cfl
      write(*,*) '    dt = ',dt
      write(*,*) '    ntime = ',ntime
      !
      !-------------------------
      ! Initialize the mesh(es)
      !-------------------------
      call init_mesh(msh(1),[-1d0,1d0],dx(1),1,order)
!      call init_mesh(msh(2),[m2start,m2start+0.75d0],dx(2),0,order)
      call initBC(msh(1)) ! only call for msh1 b/c it's the outside mesh
      !
      ! Blank out coarser overlapping cells 
      if(nmesh>1) then 
        call connect(msh(1),msh(2))
        call connect(msh(2),msh(1))
        call fixOverlap(msh(2),msh(1))
        call findIncompleteElements(msh(1),elemInfo1,nincomp1)
        call findIncompleteElements(msh(2),elemInfo2,nincomp2)

        ! cut mass matrices and merge cells if needed
        debug = 0      
        call fixMassIncompleteElements(msh(1),msh(2),elemInfo2,nincomp2,&
                consoverset,foverlap,debug)
        debug = 0
        call fixMassIncompleteElements(msh(2),msh(1),elemInfo1,nincomp1,&
                consoverset,foverlap,debug)      
      endif
      !
      ! check if cells are merged or not
      do n=1,nmesh
        isMerged(n) = 0
        do i=1,msh(n)%nelem
          if(msh(n)%parent(i).ne.i) then
            isMerged(n) = 1
            cycle
          endif
        enddo
      enddo
      write(*,*) ' '
      write(*,*) '    isMerged = ',isMerged
      !
      ! init q at t=0 
      do n=1,nmesh
        call initvar(msh(n),msh(n)%q,0)
        msh(n)%sol=msh(n)%q
      enddo
      if((nmesh>1).and.(consoverset.eq.1)) then
        call projectChild(msh(1),elemInfo1,nincomp1,msh(1)%q)
        call projectChild(msh(2),elemInfo2,nincomp2,msh(2)%q)
        msh(1)%qold = msh(1)%q 
        msh(2)%qold = msh(2)%q 
      endif
      !
      ! find exact solution f q at t=ntime*dt
      do n=1,nmesh
       call initvar(msh(n),msh(n)%qexact,ntime*dt)
      enddo
      if((nmesh>1).and.(consoverset.eq.1)) then
        call projectChild(msh(1),elemInfo1,nincomp1,msh(1)%qexact)
        call projectChild(msh(2),elemInfo2,nincomp2,msh(2)%qexact)
      endif
!      call computeMoments(msh(1),mom0(:,1),err(1),nincomp1,elemInfo1)
!      call computeMoments(msh(2),mom0(:,2),err(2),nincomp2,elemInfo2)
      do n=1,nmesh
        call output(n,msh(n))
      enddo
      write(*,*) ' '
      write(*,*) '    Initial Area Under Curve: ', sum(mom0(1,:))

      !-----------------------
      ! Iterate in time
      !-----------------------
      rk = [1d0/4d0, 8d0/15d0,5d0/12d0, 3d0/4d0];
      do i=1,ntime
!       write(*,*) '--------------------------'
!       write(*,*) 'TIMESTEP ',i
!       write(*,*) '--------------------------'
       ! RK step 1
        call timestep(nmesh,dt,msh,consoverset,elemInfo1,elemInfo2,nincomp1,nincomp2,foverlap,isupg,ireg)

        if(ieuler.eq.1) then ! Euler 1st order
          do j = 1,nmesh
            msh(j)%q=msh(j)%q+dt*msh(j)%dq
            msh(j)%sol=msh(j)%q
          enddo
          ! Project solutions onto child cells from parent cells
          if((nmesh.gt.1).and.(consoverset.eq.1)) then
             call projectChild(msh(1),elemInfo1,nincomp1,msh(1)%q)
             call projectChild(msh(2),elemInfo2,nincomp2,msh(2)%q)
          endif
          if(ilim.eq.1) then 
            if(nmesh.gt.1) then 
              write(*,*) 'Limiting Grid 1...'
              call genLimit2(msh(1)%q,msh(1)%vlim,msh(1),msh(2))
              msh(1)%q = msh(1)%vlim
              msh(1)%sol=msh(1)%q

              write(*,*) 'Limiting Grid 2...'
              call genLimit2(msh(2)%q,msh(2)%vlim,msh(2),msh(1))
              msh(2)%q = msh(2)%vlim
              msh(2)%sol=msh(2)%q
            else
              call genLimit(msh(1)%q,msh(1)%vlim,msh(1))
              msh(1)%q = msh(1)%vlim
              msh(1)%sol=msh(1)%q
            endif
          endif ! limiter
        else ! RK3
          do j = 1,nmesh
            msh(j)%q=msh(j)%sol+rk(2)*dt*msh(j)%dq
            msh(j)%sol=msh(j)%sol+rk(1)*dt*msh(j)%dq
          enddo
          ! Project solutions onto child cells from parent cells
          if((nmesh.gt.1).and.(consoverset.eq.1)) then
             call projectChild(msh(1),elemInfo1,nincomp1,msh(1)%q)
             call projectChild(msh(1),elemInfo1,nincomp1,msh(1)%sol)
             call projectChild(msh(2),elemInfo2,nincomp2,msh(2)%q)
             call projectChild(msh(2),elemInfo2,nincomp2,msh(2)%sol)
          endif
          if(ilim.eq.1) then
            if(nmesh.gt.1) then 
              call genLimit2(msh(1)%q,msh(1)%vlim,msh(1),msh(2))
              msh(1)%q = msh(1)%vlim
              call genLimit2(msh(1)%sol,msh(1)%vlim,msh(1),msh(2))
              msh(1)%sol = msh(1)%vlim
              call genLimit2(msh(2)%q,msh(2)%vlim,msh(2),msh(1))
              msh(2)%q = msh(2)%vlim
              call genLimit2(msh(2)%sol,msh(2)%vlim,msh(2),msh(1))
              msh(2)%sol = msh(2)%vlim
            else
              call genLimit(msh(1)%q,msh(1)%vlim,msh(1))
              msh(1)%q = msh(1)%vlim
              call genLimit(msh(1)%sol,msh(1)%vlim,msh(1))
              msh(1)%sol = msh(1)%vlim
            endif
          endif

          ! RK step 2
          call timestep(nmesh,dt,msh,consoverset,elemInfo1,elemInfo2,nincomp1,nincomp2,foverlap,isupg,ireg)
          do j = 1,nmesh
            msh(j)%q=msh(j)%sol+rk(3)*dt*msh(j)%dq
          enddo
          if((nmesh.gt.1).and.(consoverset.eq.1)) then
             call projectChild(msh(1),elemInfo1,nincomp1,msh(1)%q)
             call projectChild(msh(2),elemInfo2,nincomp2,msh(2)%q)
          endif
          if(ilim.eq.1) then
            if(nmesh.gt.1) then 
              call genLimit2(msh(1)%q,msh(1)%vlim,msh(1),msh(2))
              msh(1)%q = msh(1)%vlim
              call genLimit2(msh(2)%q,msh(2)%vlim,msh(2),msh(1))
              msh(2)%q = msh(2)%vlim
            else
              call genLimit(msh(1)%q,msh(1)%vlim,msh(1))
              msh(1)%q = msh(1)%vlim
            endif
          endif

          ! RK step 3
          call timestep(nmesh,dt,msh,consoverset,elemInfo1,elemInfo2,nincomp1,nincomp2,foverlap,isupg,ireg)
          do j = 1,nmesh
            msh(j)%sol=msh(j)%sol+rk(4)*dt*msh(j)%dq
          enddo
          if((nmesh.gt.1).and.(consoverset.eq.1)) then
             call projectChild(msh(1),elemInfo1,nincomp1,msh(1)%sol)
             call projectChild(msh(2),elemInfo2,nincomp2,msh(2)%sol)
          endif
          if(ilim.eq.1) then
            if(nmesh.gt.1) then 
              call genLimit2(msh(1)%sol,msh(1)%vlim,msh(1),msh(2))
              msh(1)%sol = msh(1)%vlim
              call genLimit2(msh(2)%sol,msh(2)%vlim,msh(2),msh(1))
              msh(2)%sol = msh(2)%vlim
            else
              call genLimit(msh(1)%sol,msh(1)%vlim,msh(1))
              msh(1)%sol = msh(1)%vlim
            endif
          endif

          do n = 1,nmesh
            msh(n)%q = msh(n)%sol
  
            ! copy over curr q values to qold
            msh(n)%qold = msh(n)%q
          enddo
        endif ! Euler or RK
      enddo ! timesteps
      
      ! write final output
      do n=1,nmesh
        call output(nmesh+n,msh(n))
      enddo
!      call computeMoments(msh(1),mom1(:,1),err(1),nincomp1,elemInfo1)
!      call computeMoments(msh(2),mom1(:,2),err(2),nincomp2,elemInfo2)
!      sweep(h,1) = sum(err(:))
!      sweep(h,2) = sum(mom1(1,:))
!      sweep(h,3) = sum(mom1(1,:))-sum(mom0(1,:))
!      write(*,*) '    Min Rem Frac M1: ',minval(msh(1)%dxcut)/dx(1)
!      write(*,*) '    Min Rem Frac M2: ',minval(msh(2)%dxcut)/dx(2)
      
      ! Clear memory
      do n = 1,nmesh
        call clearMem(msh(n))
      enddo
      if (allocated(elemInfo1)) deallocate(elemInfo1)
!      if (allocated(elemInfo2)) deallocate(elemInfo2)

      ! call timer
      call cpu_time(time(2))
      write(*,*) ' '
      write(*,*) '    Run time: ',time(2)-time(1)
      write(*,*) '    Run time/step: ',(time(2)-time(1))/ntime
    enddo ! mesh sweep

    write(*,*) ' '
    write(*,*) 'FINAL SWEEP OUTPUT: '
    write(*,*) '  H SWEEP L2 SQUARED ERROR = ',sweep(:,1)
    write(*,*) '  H SWEEP CONS AREA = ',sweep(:,2)
    write(*,*) '  H SWEEP CONS ERROR = ',sweep(:,3)
    write(*,*) '----------------------------------'
  enddo ! p sweep
  enddo ! overflap
  enddo ! shape function
  enddo ! cons overset
end program conservative_overset
