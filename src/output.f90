subroutine output(iout,msh)
  use code_types
  use bases
  !
  implicit none
  integer, intent(in):: iout
  type(mesh), intent(in) :: msh
  !
  integer :: i
  real*8 :: qout(msh%nshp),dqout(msh%nshp),q1,q2
  character*6 :: fname
  !
  write(fname,'(A5,I1)') 'aout.',iout
  open(unit=11,file=fname,form='formatted')
  write(fname,'(A5,I1)') 'mesh.',iout
  open(unit=12,file=fname,form='formatted')

  do i=1,msh%nelem
    if (msh%iblank(msh%e2n(1,i)) == 1 .or. &
        msh%iblank(msh%e2n(2,i)) == 1) then 
        
       call shapefunction(msh%nshp,msh%xe(1,i),[msh%xe(1,i),msh%xe(2,i)],msh%sol(1,:,i),qout,dqout)
       q1 = sum(qout)
       call shapefunction(msh%nshp,msh%xe(2,i),[msh%xe(1,i),msh%xe(2,i)],msh%sol(1,:,i),qout,dqout)
       q2 = sum(qout)

       write(11,*) 0.5*(msh%xe(1,i)+msh%xe(2,i)),&
                   (q1+q2)/msh%nshp                   
       write(12,*) msh%xe(1,i),q1
       write(12,*) msh%xe(msh%nshp,i),q2
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
