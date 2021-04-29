subroutine clearMem(msh)
  use code_types
  type(mesh), intent(inout) :: msh
  !
  if (allocated(msh%e2n)) deallocate(msh%e2n)
  if (allocated(msh%face)) deallocate(msh%face)
  if (allocated(msh%iblank)) deallocate(msh%iblank)
  if (allocated(msh%dq)) deallocate(msh%dq)
  if (allocated(msh%q)) deallocate(msh%q)
  if (allocated(msh%q0)) deallocate(msh%q0)
  if (allocated(msh%sol)) deallocate(msh%sol)
  if (allocated(msh%x)) deallocate(msh%x)
  if (allocated(msh%xe)) deallocate(msh%xe)
  if (allocated(msh%nres)) deallocate(msh%nres)
  if (allocated(msh%rhs)) deallocate(msh%rhs)
  if (allocated(msh%mass)) deallocate(msh%mass)
  if (allocated(msh%dx)) deallocate(msh%dx)
  if (allocated(msh%wgauss)) deallocate(msh%wgauss)
  if (allocated(msh%xgauss)) deallocate(msh%xgauss)
  if (allocated(msh%shp)) deallocate(msh%shp)
  if (allocated(msh%dshp)) deallocate(msh%dshp)
end subroutine clearMem
