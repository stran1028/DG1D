subroutine initvar(msh)
  !
  use code_types
  use pde
  implicit none
  !
  type(mesh), intent(inout) :: msh
  integer :: i,j
  real*8 :: xc,xx,dx
  !
  do i=1,msh%nelem
     xc=(msh%xe(1,i)+msh%xe(2,i))*0.5d0
     dx=msh%xe(2,i)-msh%xe(1,i)
     do j=1,msh%nshp
        !xx=xc+msh%xgauss(j)*dx
        xx = msh%xe(j,i)        ! Prescribing q var at node points but how to do DG?
        call initq(xx,msh%q(:,j,i))
     enddo
  enddo
  !
end subroutine initvar
        
        
  
