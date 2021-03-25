subroutine initvar(msh)
  !
  use code_types
  use pde
  use bases
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
        xx = msh%xe(j,i)    
        if(shptype.eq.'legendre') then 
          call initqleg(msh)
        else 
          call initq(xx,msh%q(:,j,i))
        endif
     enddo
  enddo
  !
end subroutine initvar
        
        
  
