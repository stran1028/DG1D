subroutine solveDQ(msh,dt)
  use code_types
  implicit none
  type(mesh), intent(inout) :: msh
  real*8, intent(in) :: dt
  integer :: i,j,k,index1,index2
  real*8 :: dx, detM, relax
  real*8, dimension(msh%nshp*msh%nshp) :: mass,L,U
  real*8, dimension(msh%nshp) :: y,b
  real*8, dimension(msh%nshp*msh%nshp) :: A
  real*8, dimension(2*msh%nshp*msh%nshp) :: Areg,AregT
  real*8, dimension(2*msh%nshp) :: breg
  real*8:: lambda,tmp1(8),tmp2(12),tmp3(6) 
  !
  do i=1,msh%nelem
    ! solve Ax=b problem for x

    ! use Tikhonov Regularization on cut elements for stability
    if(msh%dxcut(i)/msh%dx(i).lt.0.50) then ! only do for severely cut elements 
      lambda = 0d0    ! Adding the least amount of bias as possible

      ! assemble the regularized matrix [A; lambdaI] and vector [b; 0]
      breg(1:msh%nshp) = msh%rhs(1,:,i)
      breg(msh%nshp+1:2*msh%nshp) = 0d0 

      Areg = 0d0
      AregT = 0d0
      do j = 1,2*msh%nshp
      do k = 1,msh%nshp
         index1 = (j-1)*msh%nshp + k 
         index2 = (k-1)*2*msh%nshp + j ! Note j and k are switched b/c I want the transpose
         if(index1.le.msh%nshp*msh%nshp) then 
           Areg(index1) = msh%mass(1,index1,i)
         else
           if(j-msh%nshp.eq.k) Areg(index1) = lambda
           if(j-msh%nshp.ne.k) Areg(index1) = 0d0
         endif
         AregT(index2) = Areg(index1)
      enddo
      enddo

      ! compute trans(A)*A and trans(A)*b
      call matmat(AregT,Areg,A,msh%nshp,2*msh%nshp,msh%nshp)
      call matvec(AregT,breg,msh%nshp,2*msh%nshp,b)

      !write(*,*) ' '
      !write(*,*) 'eid ',i,msh%dxcut(i)/msh%dx(i),msh%dxcut(i),msh%dx(i)
      !write(*,*) 'x ',msh%x(msh%e2n(:,i))
      !write(*,*) 'iblank ',msh%iblank(:,i)
      !write(*,*) 'mass0 ',msh%mass(1,:,i)
      !write(*,*) 'rhs0 ',msh%rhs(1,:,i)
      !write(*,*) 'mass1 ',A
      !write(*,*) 'rhs1 ',b
      !write(*,*) ' '

    else ! otherwise, do nothing
      A = msh%mass(1,:,i)
      b = msh%rhs(1,:,i)

    endif

    ! Solve Ax=b using LU decomp
    call lu(A,msh%nshp,L,U)
    call forwprop(L,b,msh%nshp,y)
    call backprop(U,y,msh%nshp,msh%dq(1,:,i))

  enddo
  !
end subroutine solveDQ
