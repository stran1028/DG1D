  subroutine matmat(A,B,C,p,q,r)
    implicit none
    integer, intent(in) :: p,q,r
    real*8, dimension(p*q),intent(in) ::  A
    real*8, dimension(q*r),intent(in) ::  B
    real*8, dimension(p*r),intent(inout) :: C
    real*8 :: tmp
    integer i,j,k,aa,bb,cc

    ! Matrix A is p x q (rows x col)
    ! Matrix B is q x r 
    do i = 1,p ! rows of C
    do j = 1,r ! cols of C
      cc = (i-1)*r+j
      tmp = 0d0
      do k = 1,q ! cols of A and rows of B
        aa = (i-1)*q+k
        bb = (k-1)*r+j
        tmp = tmp + A(aa)*B(bb) 
      enddo
      C(cc) = tmp
    enddo 
    enddo 
 
  end subroutine matmat

  subroutine matvec(A,x,m,n,b)
    implicit none
    integer, intent(in) :: m,n
    real*8, dimension(m*n),intent(in) ::  A
    real*8, dimension(n),intent(in) ::  x
    real*8, dimension(m),intent(inout) ::  b
    integer :: i,j,k

    ! m = rows in A
    ! n = cols in A, rows, in x
    b = 0d0
    do i = 1,m
      do j = 1,n
        b(i) = b(i) + A((i-1)*n+j)*x(j)
      enddo ! j loop
    enddo ! i loop

  end subroutine matvec

  subroutine backprop(U,b,n,x)
    implicit none
    integer, intent(in) :: n
    real*8, dimension(n*n),intent(in) ::  U
    real*8, dimension(n),intent(in) ::  b
    real*8, dimension(n),intent(inout) ::  x
    integer :: i,j,k,index1
  
    ! solve Ux = b where U is upper triangular matrix
    x(n) = b(n)/(U(n*n)+1e-15)
    do i = n-1,1,-1
      x(i) = b(i)
      do j = n,i+1,-1
        index1 = (i-1)*n+j
        x(i) = x(i)-x(j)*U(index1)
      enddo
      index1 = (i-1)*n+i
      x(i) = x(i)/(U(index1)+1e-15)
    enddo
  end subroutine backprop

  subroutine forwprop(L,b,n,x)
    implicit none
    integer, intent(in) :: n
    real*8, dimension(n*n),intent(in) ::  L
    real*8, dimension(n),intent(in) ::  b
    real*8, dimension(n),intent(inout) ::  x
    integer :: i,j,k,index1
  
    ! solve Lx = b where L is lower triangular matrix
    x(1) = b(1)/(L(1)+1e-15)

    do i = 2,n
      x(i) = b(i)
      do j = 1,i-1
        index1 = (i-1)*n+j
        x(i) = x(i)-x(j)*L(index1)
      enddo
      index1 = (i-1)*n+i
      x(i) = x(i)/(L(index1)+1e-15)
    enddo
  end subroutine forwprop

  subroutine lu(A,n,L,U)
    implicit none
    integer, intent(in) :: n
    real*8, dimension(n*n),intent(in) ::  A
    real*8, dimension(n*n),intent(inout) ::  L,U
    integer :: i,j,k,index1,index2
    real*8 :: s

    L = 0d0
    u = 0d0
    do i = 1,n
        ! upper matrix
        do k = i,n
          s = 0d0
          do j = 1,i
            index1 = (i-1)*n+j
            index2 = (j-1)*n+k
            s = s + l(index1)*u(index2)
          enddo ! loop j
          index1 = (i-1)*n+k
          u(index1) = A(index1)-s
        enddo ! loop k

        ! Lower matrix
        do k = i,n
          if(i.eq.k) then
            L((i-1)*n+i) = 1d0
          else
            s = 0d0
            do j = 1,i
              index1 = (k-1)*n+j
              index2 = (j-1)*n+i
              s = s + l(index1)*u(index2)
            enddo !loop j
            index1 = (k-1)*n+i
            index2 = (i-1)*n+i
            l(index1) = (A(index1)-s)/u(index2)
          endif
        enddo ! loop k
    enddo ! loop i

  end subroutine lu
