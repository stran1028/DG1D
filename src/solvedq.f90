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
    call lu(msh%mass(1,:,i),msh%nshp,L,U)
    call forwprop(L,msh%rhs(1,:,i),msh%nshp,y)
    call backprop(U,y,msh%nshp,msh%dq(1,:,i))

  enddo
  !
end subroutine solveDQ
