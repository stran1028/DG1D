subroutine timestep(nmesh,dt,msh)
   use code_types
   use pde
   use overset
   implicit none

   real*8, allocatable :: elemInfo1(:),elemInfo2(:)
   integer, intent(in) :: nmesh
   real*8, intent(in) :: dt
   type(mesh), intent(inout) :: msh(nmesh)
   integer :: n,nincomp1,nincomp2

   do n=1,nmesh
      call computeRHS(msh(n))
   enddo
   ! Adjust RHS based on overlaps
   if (nmesh > 1) then
      call findIncompleteElements(msh(1),elemInfo1,nincomp1)
      call findIncompleteElements(msh(2),elemInfo2,nincomp2)
write(*,*) 'mesh 1 to 2'
      call fixfluxIncompleteElements(msh(1),msh(2),elemInfo2,nincomp2)
write(*,*) 'mesh 2 to 1'
      call fixfluxIncompleteElements(msh(2),msh(1),elemInfo1,nincomp1)
   end if
   do n=1,nmesh
      ! Solve M dq = rhs
      call solvedq(msh(n),dt)
!      call computeMoments(msh(n),mom1(:,n))
   enddo
end subroutine timestep
