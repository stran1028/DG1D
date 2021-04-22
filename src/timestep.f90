subroutine timestep(nmesh,dt,msh,consoverset,elemInfo1,elemInfo2,nincomp1,nincomp2)
   use code_types
   use pde
   use overset
   implicit none

   integer, intent(in) :: nmesh,consoverset
   integer :: n,nincomp1,nincomp2
   real*8, intent(inout) :: elemInfo1(nincomp1),elemInfo2(nincomp2)
   real*8, intent(in) :: dt
   type(mesh), intent(inout) :: msh(nmesh)
   do n=1,nmesh
      call computeRHS(msh(n))
   enddo

   ! Adjust RHS based on overlaps
   if (nmesh > 1) then
      call fixfluxIncompleteElements(msh(1),msh(2),elemInfo2,nincomp2,consoverset)
      call fixfluxIncompleteElements(msh(2),msh(1),elemInfo1,nincomp1,consoverset)
   end if
!   write(*,*) 'debug: ref mass  = ',msh(2)%mass(1,:,3)
!   write(*,*) 'debug: mass e1 = ',msh(2)%mass(1,:,1)
!   write(*,*) 'debug: rhs e1 = ',msh(2)%rhs(1,:,1)
   do n=1,nmesh
      ! Solve M dq = rhs
      call solvedq(msh(n),dt)
   enddo
end subroutine timestep
