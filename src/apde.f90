module pde
  character*20 :: pde_descriptor
  integer, save :: pde_set=0
  real*8 :: a
contains
  subroutine set_type(pdetype,wavespeed)
    implicit none
    character*(*) :: pdetype
    real*8, optional, intent(in) :: wavespeed
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
  end subroutine set_type
  !
  subroutine initqleg(msh)
    use code_types
    use bases
    type(mesh), intent(inout)::msh
    integer i,j,k
!    real*8 :: mass(4),detM,invM(4)
    real*8,dimension(msh%nshp*msh%nshp) :: L,U
    real*8,dimension(msh%nshp) :: y
    real*8 :: qvals(msh%nshp),dqvals(msh%nshp),f(msh%nshp),tmp,xloc,dx

    ! Use projection to set the IC
    ! int(NaNbub) = int(Na y(x,0))
    ! ub = M^-1 int(Na y(x,0))
    do i = 1,msh%nelem
      dx = msh%xe(2,i)-msh%xe(1,i)

!      mass = [msh%mass(1,1,1,i), msh%mass(1,1,2,i), msh%mass(1,2,1,i),msh%mass(1,2,2,i)]
!      detM = mass(1)*mass(4) - mass(2)*mass(3)
!      invM = [mass(4),-mass(2),-mass(3),mass(1)]
!      invM = invM/(1e-16+detM)

      f = 0.0d0
      do j = 1,msh%ngauss
        xloc = (msh%xgauss(j)+0.5d0)*dx + msh%xe(1,i)
        tmp = exp(-0.01d0*xloc*xloc)
        if (abs(xloc) > 40d0) tmp=0d0
!        tmp = exp(-20d0*xloc*xloc)
!        if (abs(xloc) > 0.8d0) tmp=0d0
        !tmp = exp(-.010d0*xloc*xloc)
        !if (abs(xloc) > 80d0) tmp=0d0
        call shapefunction(msh%nshp,msh%xgauss(j),[-0.5d0,0.5d0],[1d0,1d0],qvals,dqvals)
        do k = 1,msh%nshp
          f(k) = f(k) + msh%wgauss(j)*msh%dx(i)*tmp*qvals(k)
        enddo
      enddo

      ! solve the elementary system for u
      call lu(msh%mass,msh%nshp,L,U)
      call backpropL(L,f,msh%nshp,y)
      call backpropU(U,y,msh%nshp,msh%q(1,:,i))

!      msh%q(1,1,i) = invM(1)*f(1) + invM(2)*f(2)
!      msh%q(1,2,i) = invM(3)*f(1) + invM(4)*f(2)
     
    enddo
  end subroutine initqleg
  !
  subroutine initq(x,q)
    implicit none
    real*8, intent(in) :: x
    real*8, intent(out) :: q(:)
!    q=1d0
    !return

    if ((index(pde_descriptor,'linear_advection') .le. 0) .and. &
        (index(pde_descriptor, 'burgers') .le. 0)) then
       write(6,*) 'problem type not implemented'
       call exit(1)
    else     
!       q = 1d0
!        q = sin(x)
!       q=exp(-20*x*x)
!       if (abs(x) > 0.8d0) q=0d0
       q=exp(-0.01*x*x)
       if (abs(x) > 40d0) q=0d0
    endif
  end subroutine initq
  !
  subroutine flux(ql,qr,flx)
    !
    implicit none
    !
    real*8, intent(in) :: ql
    real*8, intent(in) :: qr
    real*8, intent(out) :: flx
    !

    ! Lax Friedrichs Flux
    if (index(pde_descriptor,'linear_advection') > 0 ) then
     flx=0.5d0*(a*(ql+qr)-abs(a)*(qr-ql))
    else if (index(pde_descriptor,'burgers') > 0) then
     flx=0.5d0*(0.5d0*(ql*ql+qr*qr)-0.5d0*(ql+qr)*(qr-ql))
    endif
  end subroutine flux
  !
  subroutine volintstrong(dqtmp,vol)
    !
    implicit none
    !
    real*8, intent(in) :: dqtmp
    real*8, intent(out) :: vol
    !
    if (index(pde_descriptor,'linear_advection') > 0 ) then
     vol=a*dqtmp
    else if (index(pde_descriptor,'burgers') > 0) then
     write(*,*) 'Whoops. I will add this later'
     call exit(1)
    endif
  end subroutine volintstrong
  !
  subroutine volint(qtmp,vol)
    !
    implicit none
    !
    real*8, intent(in) :: qtmp
    real*8, intent(out) :: vol
    !
    if (index(pde_descriptor,'linear_advection') > 0 ) then
     vol=a*qtmp
    else if (index(pde_descriptor,'burgers') > 0) then
     vol=0.5*qtmp*qtmp
    endif
  end subroutine volint
  !
end module pde
