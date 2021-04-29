subroutine output(iout,msh)
  use code_types
  use pde
  use bases
  !
  implicit none
  integer, intent(in):: iout
  type(mesh), intent(in) :: msh
  !
  integer :: i,j
  real*8 :: qout(msh%nshp),dqout(msh%nshp),q1,q2,error(1)
  character*22 :: fname
  !
  write(fname,'(A7,I0.3)') 'iblank.',iout
  open(unit=11,file=fname,form='formatted')
  write(fname,'(A5,I0.3)') 'mesh.',iout
  open(unit=12,file=fname,form='formatted')

  do i=1,msh%nelem
    write(11,*) i,msh%xe(1,i),msh%xe(2,i),msh%iblank(:,i)

    if (maxval(msh%iblank(:,i)) == 1) then 
       do j = 1,msh%nshp        
         call shapefunction(msh%nshp,msh%x(msh%e2n(j,i)),[msh%xe(1,i),msh%xe(2,i)],msh%sol(1,:,i),qout,dqout)
         q1 = sum(qout)

         call initq(msh%x(msh%e2n(j,i)),error)
         error = error-q1
         write(12,*) msh%x(msh%e2n(j,i)),q1, error
      enddo

    endif
  enddo
  close(11)
  close(12)
  !
end subroutine output
