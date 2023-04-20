module pde
  character*20 :: pde_descriptor
  integer, save :: pde_set=0
  real*8 :: a,diff
contains
  subroutine set_type(pdetype,wavespeed,muinf)
    implicit none
    character*(*) :: pdetype
    real*8, optional, intent(in) :: wavespeed,mu
    write(pde_descriptor,'(A20)') pdetype
    pde_descriptor=trim(adjustl(pde_descriptor))
    if (pde_descriptor .ne. 'linear_advection' .and. & 
        pde_descriptor .ne. 'burgers' ) then
       write(6,*) 'Only linear advection and burgers are implemented'
       stop
    endif
    pde_set=1
    if (present(wavespeed)) then
       a=wavespeed
    else
       a=1d0
    endif
    if (present(muinf)) then
       mu=muinf
    else
       mu=0d0
    endif

  end subroutine set_type
  !
  subroutine initqleg(msh)
    use code_types
    use bases
    type(mesh), intent(inout)::msh
    integer i,j,k
    real*8,dimension(msh%nshp*msh%nshp) :: L,U
    real*8,dimension(msh%nshp) :: y
    real*8 :: qvals(msh%nshp),dqvals(msh%nshp),f(msh%nshp),tmp,xloc,dx

    ! Use projection to set the IC
    ! int(NaNbub) = int(Na y(x,0))
    ! ub = M^-1 int(Na y(x,0))
    do i = 1,msh%nelem
      dx = msh%xe(2,i)-msh%xe(1,i)
      f = 0.0d0
      do j = 1,msh%ngauss
        xloc = (msh%xgauss(j)+0.5d0)*dx + msh%xe(1,i)
!        tmp = exp(-0.01d0*xloc*xloc)
!        if (abs(xloc) > 40d0) tmp=0d0
        tmp = 0.1d0*exp(-20d0*xloc*xloc)
        if (abs(xloc) > 0.8d0) tmp=0d0
        !tmp = exp(-.010d0*xloc*xloc)
        !if (abs(xloc) > 80d0) tmp=0d0
        qvals = 1d0
        call shapefunction(msh%nshp,msh%xgauss(j),[-0.5d0,0.5d0],qvals,qvals,dqvals)
        do k = 1,msh%nshp
          f(k) = f(k) + msh%wgauss(j)*msh%dx(i)*tmp*qvals(k)
        enddo
      enddo

      ! solve the elementary system for u
      call lu(msh%mass(1,:,i),msh%nshp,L,U)
      call forwprop(L,f,msh%nshp,y)
      call backprop(U,y,msh%nshp,msh%q(1,:,i))

    enddo
  end subroutine initqleg
  !
  subroutine initq(x,q)
    implicit none
    real*8, intent(in) :: x
    real*8, intent(out) :: q(:)

    if ((index(pde_descriptor,'linear_advection') .le. 0) .and. &
        (index(pde_descriptor, 'burgers') .le. 0)) then
       write(6,*) 'problem type not implemented'
       call exit(1)
    else     
!       q = 1d0
!       if (x < 1d0) q=0d0
!       if (x > 5d0) q=0d0


       q=exp(-20*x*x)
       if (abs(x) > 0.8d0) q=0d0
!       q=exp(-0.01*x*x)
!       if (abs(x) > 40d0) q=0d0
    endif
  end subroutine initq
  !
  subroutine flux(ql,qr,dql,dqr,flx)
    !
    implicit none
    !
    real*8, intent(in) :: ql,dql
    real*8, intent(in) :: qr,dqr
    real*8, intent(out) :: flxa,flxd
    real*8 :: C11,C12
    !

    if (index(pde_descriptor,'linear_advection') > 0 ) then
      ! Lax Friedrichs Flux
      ! flxa=0.5d0*(a*(ql+qr)-abs(a)*(qr-ql))

      C11 = 0.0d0
      C12 = 0.5d0 ! 0.5 dot n-

      flx=a*(0.5d0*(ql+qr) + C12*(ql-qr) )
      dflx = mu*(0.5d0*(dql+dqr) + C11*(ql-qr) - C12*(dql-dqr))

      flx=flxa-flxd
    else if (index(pde_descriptor,'burgers') > 0) then
     flx=0.5d0*(0.5d0*(ql*ql+qr*qr)-0.5d0*(ql+qr)*(qr-ql))
    endif
  end subroutine flux
  !
  subroutine volint(qtmp,dqtmp,vol)
    !
    implicit none
    !
    real*8, intent(in) :: qtmp,dqtmp
    real*8, intent(out) :: vol
    !
    if (index(pde_descriptor,'linear_advection') > 0 ) then
     vol=a*qtmp - mu*dqtmp
    else if (index(pde_descriptor,'burgers') > 0) then
     vol=0.5*qtmp*qtmp
    endif
  end subroutine volint
  !
end module pde
