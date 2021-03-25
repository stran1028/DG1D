subroutine timeIntegrate(msh,dt)
  use code_types
  implicit none
  type(mesh), intent(inout) :: msh
  real*8, intent(in) :: dt
  integer :: i,j
  real*8 :: dx, detM, relax
  real*8, dimension(2) :: dudt
  real*8, dimension(4) :: inVM
  !
  do i=1,msh%nelem
    dx=msh%xe(2,i)-msh%xe(1,i)

    ! Hardcoding for nshp=2 for linear advection for now
    
    ! Need to multiply mass matrix by dx since it wasn't done yet
    detM = (msh%mass(1,1,1,i)*msh%mass(1,2,2,i) - msh%mass(1,1,2,i)*msh%mass(1,2,1,i))
    invM = [msh%mass(1,2,2,i),-msh%mass(1,1,2,i),-msh%mass(1,2,1,i),msh%mass(1,1,1,i)]
    invM = invM/(1e-16+detM)

    dudt(1) = invM(1)*msh%rhs(1,1,i) + invM(2)*msh%rhs(1,2,i)
    dudt(2) = invM(3)*msh%rhs(1,1,i) + invM(4)*msh%rhs(1,2,i)

    if(i.eq.50) then 
    write(*,*) 'Element ',i
    write(*,*) '  mass: ',msh%mass(1,1,1,i),msh%mass(1,1,2,i),msh%mass(1,2,1,i),msh%mass(1,2,2,i)
    write(*,*) '  detM: ',detM
    write(*,*) '  invmass: ',invM(1),invM(2),invM(3),invM(4)
    write(*,*) '  rhs: ',msh%rhs(1,:,i)
    write(*,*) '  dudt: ',dudt(1),dudt(2)
    write(*,*) '  u0: ',msh%q(1,1,i),msh%q(1,2,i)
    endif
   ! 1st order Euler in time
    relax = 1.0
    msh%q(1,1,i) = msh%q(1,1,i) + relax*dudt(1)*dt
    msh%q(1,2,i) = msh%q(1,2,i) + relax*dudt(2)*dt
    if(i.eq.50) write(*,*) '  unew: ',msh%q(1,1,i),msh%q(1,2,i)
  enddo
  !
end subroutine timeIntegrate
