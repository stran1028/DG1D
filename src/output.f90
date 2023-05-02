subroutine output(iout,msh)
  use code_types
  use pde
  use bases
  !
  implicit none
  integer, intent(in):: iout
  type(mesh), intent(in) :: msh
  !
  integer :: i,j,npts,bound
  real*8 :: qtmp(msh%nshp),dqtmp(msh%nshp),q1,ic,error(1),xloc,dx
  character*22 :: fname
  !
  write(fname,'(A7,I0.3)') 'iblank.',iout
  open(unit=11,file=fname,form='formatted')
  write(fname,'(A5,I0.3)') 'mesh.',iout
  open(unit=12,file=fname,form='formatted')

  npts = 2*msh%nshp
  do i=1,msh%nelem
    write(11,*) i,msh%xe(1,i),msh%xe(2,i),msh%iblank(:,i)

    if (maxval(msh%iblank(:,i)) == 1) then 
!       do j = 1,msh%nshp       
       ! output data at npts equidistant points along element
       do j = 1,npts 
         if((j.eq.1).or.(j.eq.npts)) then 
           bound = 1
         else
           bound = 0
         endif
         dx = msh%dx(i)/(npts-1)
         xloc = msh%xe(1,i) + (j-1)*dx
         q1 = sum(qtmp)

         call shapefunction(msh%nshp,xloc,[msh%xe(1,i),msh%xe(2,i)],msh%sol(1,:,i),qtmp,dqtmp)
         q1 = sum(qtmp)
         call shapefunction(msh%nshp,xloc,[msh%xe(1,i),msh%xe(2,i)],msh%qexact(1,:,i),qtmp,dqtmp)
         ic = sum(qtmp)
         error = q1-ic
         write(12,*) xloc,q1, error,bound

!         call shapefunction(msh%nshp,msh%x(msh%e2n(j,i)),[msh%xe(1,i),msh%xe(2,i)],msh%sol(1,:,i),qtmp,dqtmp)
!         q1 = sum(qtmp)
!         error = msh%qexact(1,j,i)-q1
!         write(12,*) msh%x(msh%e2n(j,i)),q1, error
      enddo

    endif
  enddo
  close(11)
  close(12)
  !
end subroutine output
