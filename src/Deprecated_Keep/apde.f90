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
  subroutine initq(x,q)
    implicit none
    real*8, intent(in) :: x
    real*8, intent(out) :: q(:)
    !q=1d0
    !return
    q = 0d0
    if (index(pde_descriptor,'linear_advection') > 0 .or. &
        index(pde_descriptor, 'burgers') > 0) then
        q = exp(-0.001*x*x)
!       q=exp(-20*x*x)
!       if (abs(x) > 0.8d0) q=0d0
       
!       q = 1-abs(x)/50d0
!       q = maxval([q,0d0])
      
    else
       write(6,*) 'not implemented'
       stop
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
