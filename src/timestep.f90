subroutine timestep(nmesh,dt,msh,consoverset,elemInfo1,elemInfo2,nincomp1,nincomp2,foverlap,isupg,ireg)
   use code_types
   use pde
   use overset
   implicit none

   integer, intent(in) :: nmesh,consoverset,isupg,ireg
   integer :: n,nincomp1,nincomp2,debug
   real*8, intent(inout) :: elemInfo1(3,nincomp1),elemInfo2(3,nincomp2)
   real*8, intent(inout) :: foverlap
   real*8, intent(in) :: dt
   type(mesh), intent(inout) :: msh(nmesh)

   ! compute interior element RHS
   do n=1,nmesh
      call computeRHS(msh(n),isupg,dt)
   enddo

   ! Adjust RHS based on overlaps
   if (nmesh > 1) then
      debug = 0
      call fixfluxIncompleteElements(msh(1),msh(2),elemInfo2,nincomp2,&
              consoverset,foverlap,isupg,dt,debug)
      debug = 0
      call fixfluxIncompleteElements(msh(2),msh(1),elemInfo1,nincomp1,&
              consoverset,foverlap,isupg,dt,debug)
   end if

   do n=1,nmesh
      ! Solve M dq = rhs
      call solvedq(msh(n),dt,ireg)
   enddo

!   ! Project solutions onto child cells from parent cells
!   if((nmesh.gt.1).and.(consoverset.eq.1)) then 
!      call projectChild(msh(1),elemInfo1,nincomp1)
!      call projectChild(msh(2),elemInfo2,nincomp2)
!   endif

end subroutine timestep
