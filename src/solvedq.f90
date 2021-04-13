subroutine solveDQ(msh,dt)
  use code_types
  implicit none
  type(mesh), intent(inout) :: msh
  real*8, intent(in) :: dt
  integer :: i,j
  real*8 :: dx, detM, relax
  real*8, dimension(4) :: inVM,mass
  !
  do i=1,msh%nelem
    ! Hardcoding for nshp=2 for linear advection for now
    
    mass = [msh%mass(1,1,1,i), msh%mass(1,1,2,i), msh%mass(1,2,1,i),msh%mass(1,2,2,i)]
    detM = mass(1)*mass(4) - mass(2)*mass(3)
    invM = [mass(4),-mass(2),-mass(3),mass(1)]
    invM = invM/(1e-16+detM)
    !msh%rhs(:,:,i) = -msh%rhs(:,:,i)
    msh%dq(1,1,i) = invM(1)*msh%rhs(1,1,i) + invM(2)*msh%rhs(1,2,i)
    msh%dq(1,2,i) = invM(3)*msh%rhs(1,1,i) + invM(4)*msh%rhs(1,2,i)
 
    if(i.eq.131) then 
    write(*,*) 'Element ',i
    write(*,*) '  mass = ',mass(1),mass(2),mass(3),mass(4)
    write(*,*) '  iblank = ',msh%iblank(msh%e2n(:,i))
    write(*,*) '  x = ',msh%x(msh%e2n(:,i))
    write(*,*) '  dx = ',msh%dx(i)
    write(*,*) '  detM = ',detM
    write(*,*) '  invM = ',invM
    write(*,*) '  rhsV = ',msh%rhsv(:,:,i)   
    write(*,*) '  rhsF = ',msh%rhsf(:,:,i)   
    write(*,*) '  rhs = ',msh%rhs(:,:,i)   
    write(*,*) '  q0 = ',msh%q(1,:,i) 
    write(*,*) '  msh%dq = ',msh%dq(1,:,i)
    endif

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
