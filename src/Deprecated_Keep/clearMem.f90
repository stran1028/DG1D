subroutine clearMem(msh)
  use code_types
  type(mesh), allocatable :: msh
  !
  if (allocated(msh%e2n)) deallocate(msh%e2n)
  if (allocated(msh%face)) deallocate(msh%face)
  if (allocated(msh%iblank)) deallocate(msh%iblank)
  if (allocated(msh%q)) deallocate(msh%q)
  if (allocated(msh%x)) deallocate(msh%x)
  if (allocated(msh%xe)) deallocate(msh%xe)
  if (allocated(msh%nres)) deallocate(msh%nres)
  if (allocated(msh%rhs)) deallocate(msh%rhs)
end subroutine clearMem
