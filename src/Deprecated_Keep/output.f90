subroutine output(iout,msh)
  use code_types
  !
  implicit none
  integer, intent(in):: iout
  type(mesh), intent(in) :: msh
  !
  integer :: i
  character*6 :: fname
  !
  write(fname,'(A5,I1)') 'aout.',iout
  open(unit=11,file=fname,form='formatted')
  do i=1,msh%nelem
    if (msh%iblank(msh%e2n(1,i)) == 1 .or. &
        msh%iblank(msh%e2n(2,i)) == 1) then 
       write(11,*) 0.5*(msh%xe(1,i)+msh%xe(2,i)),&
                   0.5*(msh%q(1,1,i)+msh%q(1,2,i))
    endif
  enddo
  close(11)
  !
  write(fname,'(A5,I1)') 'mesh.',iout
  open(unit=11,file=fname,form='formatted')
  do i=1,msh%nelem
    if (msh%iblank(msh%e2n(1,i)) == 1 .or. &
        msh%iblank(msh%e2n(2,i)) == 1) then 
       write(11,*) msh%xe(1,i),msh%q(1,1,i)
       write(11,*) msh%xe(2,i),msh%q(1,2,i)
    endif
  enddo
  close(11)
!  do i=1,msh%nnodes
!    write(20+iout,*) msh%iblank(i)
!  enddo
  !
end subroutine output
