subroutine initvar(msh,vec,t)
  !
  use code_types
  use pde
  use bases
  implicit none
  !
  type(mesh), intent(inout) :: msh
  real*8, intent(inout) :: vec(msh%nfields,msh%nshp,msh%nelem)
  real*8, intent(in) :: t
  integer :: i,j
  real*8 :: xx
  real*8 :: tmp(msh%nfields,msh%nshp,msh%nelem)
  !
  do i=1,msh%nelem
     do j=1,msh%nshp
        if(shptype.eq.'legendre') then 
          call initqleg(msh,tmp,t)
        else 
          xx = msh%x(msh%e2n(j,i))
          call initq(xx,tmp(:,j,i),t)
        endif
     enddo
  enddo
  !
  ! Copy over to output
  vec = tmp
  !
end subroutine initvar
