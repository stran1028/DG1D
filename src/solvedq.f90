subroutine solveDQ(msh,dt)
  use code_types
  implicit none
  type(mesh), intent(inout) :: msh
  real*8, intent(in) :: dt
  integer :: i,j
  real*8 :: dx, detM, relax
  real*8, dimension(msh%nshp*msh%nshp) :: mass,L,U
  real*8, dimension(msh%nshp) :: y
  !
  do i=1,msh%nelem
    ! solve Ax=b problem for x
    call lu(msh%mass,msh%nshp,L,U)
    call backpropL(L,msh%rhs(1,:,i),msh%nshp,y)
    call backpropU(U,y,msh%nshp,msh%dq(1,:,i))

    write(*,*) ' '
    write(*,*) 'i = ',i
    write(*,*) 'mass = ',msh%mass(1,:,i)
    write(*,*) 'rhs = ',msh%rhs(1,:,i)
    write(*,*) 'dq = ',msh%dq(1,:,i)
    
   ! 1st order Euler in time
!    relax = 1.0d0
!    msh%q(1,1,i) = msh%q(1,1,i) + relax*msh%dq(1,i)*dt
!    msh%q(1,2,i) = msh%q(1,2,i) + relax*msh%dq(2,i)*dt
    
!    if(abs(i-50).lt.10) then 
!    write(*,*) '  q1 = ',msh%q(:,:,i) 
!    write(*,*) ' '
!    endif

  enddo
  !
end subroutine solveDQ
