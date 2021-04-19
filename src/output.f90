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
  character*6 :: fname
  !
  write(fname,'(A5,I1)') 'iblank.',iout
  open(unit=11,file=fname,form='formatted')
  write(fname,'(A5,I1)') 'mesh.',iout
  open(unit=12,file=fname,form='formatted')

  do i=1,msh%nelem
    write(11,*) msh%xe(1,i),msh%xe(2,i),msh%iblank(:,i)
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
!  do i=1,msh%nnodes
!    write(20+iout,*) msh%iblank(i)
!  enddo
  !
end subroutine output
