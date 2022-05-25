module slopelimiter
contains
  ! modified minmod function (Warburton 153)
  subroutine minmod2(v,vout,m,dx)
    implicit none
    integer, intent(in) :: m    ! number of arguments into min mod
    real*8, intent(in)  :: v(m),dx
    real*8, intent(inout) :: vout

    integer :: i,j,k
    real*8 :: fact ! tuning factor, low is higher dissipation
    real*8 :: a(m)
    
    fact = 100d0*dx*dx
!    a(1) = v(1)
!    do i = 2,m
!      a(i) = v(i) + fact*dx*dx*(abs(v(i))/v(i))
!    enddo
    if(abs(v(1)).gt.fact) then 
      call minmod(v,vout,m)
    else
      vout = v(1)
    endif

  end subroutine

  ! minmod function (see Warburton pg 150)
  subroutine minmod(v,vout,m)
    
    implicit none
    integer, intent(in) :: m    ! number of arguments into min mod
    real*8, intent(in)  :: v(m)
    real*8, intent(inout) :: vout

    integer :: i,j,k
    real*8 :: s
    
    s = 0
    do i = 1,m
      s = s + v(i)/abs(v(i))
    enddo
    s = s/m

    if(abs(s).eq.1) then ! if slopes all same sign, return min. mag slope
      vout = s*minval(abs(v))
    else      ! if slopes switch signs, set slope to 0
      vout = 0.0d0; 
    endif
  end subroutine minmod
  !
  ! Pi^k Slope limiter from Cockburn (or Warburton pg 151) 
  subroutine genLimit2(vec,vecout,mshA,mshB)
    use code_types
    use pde
    use bases
    implicit none
    !
    type(mesh), intent(in) :: mshA,mshB
    real*8,intent(in) :: vec(mshA%nfields,mshA%nshp,mshA%nelem)
    real*8,intent(inout) :: vecout(mshA%nfields,mshA%nshp,mshA%nelem)
    real*8 :: u0, u0m1, u0p1
    real*8 :: xg, xb1, xb2
    real*8 :: du0, du0m1, du0p1
    real*8 :: uL, uR, uL2, uR2, du
    real*8 :: qvals(mshA%nshp), dqvals(mshA%nshp)
    real*8 :: dh,tmp,f(mshA%nshp),diff,tol,dx
    real*8,dimension(mshA%nshp*mshA%nshp) :: mass,L,U
    real*8,dimension(mshA%nshp) :: y
    integer :: n,i,j,k,eL,eR,idebug,aa,bb,index1
    do i = 1,mshA%nelem
      eL=mshA%face(1,i) ! L neigh elem
      eR=mshA%face(2,i) ! R neigh elem
      dx = mshA%dx(i)
    idebug = 0
      do n = 1,mshA%nfields
        if(maxval(mshA%iblank(:,i)).eq.1) then ! ignore blanked cells
          if(minval(mshA%iblank(:,i)).eq.1) then ! interior element
            ! Convert elements to linear and 
            ! get the cell averaged quantities u0 for elements i-1, i, and i+1
            u0m1 = 0.0d0
            u0p1 = 0.0d0
!            if(shptype.eq.'legendre') then 
              do j = 1,mshA%ngauss
                ! need to get from other mesh if on boundary
                ! how does this work for overset?
                qvals = 0d0
                dqvals = 0d0
                call shapefunction(2,mshA%xgauss(j),[-0.5d0,0.5d0],vec(n,:,eL),qvals,dqvals)
                u0m1 = u0m1 + SUM(qvals)*mshA%wgauss(j) !jacobian and 1/area cancel out
                call shapefunction(2,mshA%xgauss(j),[-0.5d0,0.5d0],vec(n,:,eR),qvals,dqvals)
                u0p1 = u0p1 + SUM(qvals)*mshA%wgauss(j)
              enddo

!            elseif(shptype.eq.'lagrange') then 
!              do j = 1,mshA%ngauss
!                call shapefunction(mshA%nshp,mshA%xgauss(j),[-0.5d0,0.5d0],vec(n,:,eL),qvals,dqvals)
!                u0m1 = u0m1 + SUM(qvals)*mshA%wgauss(j) !jacobian and 1/area cancel out
!                call shapefunction(mshA%nshp,mshA%xgauss(j),[-0.5d0,0.5d0],vec(n,:,i),qvals,dqvals)
!                u0   = u0   + SUM(qvals)*mshA%wgauss(j)
!                du0   = du0 + SUM(dqvals)*mshA%wgauss(j)
!                call shapefunction(mshA%nshp,mshA%xgauss(j),[-0.5d0,0.5d0],vec(n,:,eR),qvals,dqvals)
!                u0p1 = u0p1 + SUM(qvals)*mshA%wgauss(j)
!              enddo
!
!              ! approximate un-limited end points 
!              uL = u0-du0*dx*0.50d0
!              uR = u0+du0*dx*0.50d0
!            endif ! shape functions
          else if(mshA%iblank(1,i).ne.1) then ! leftmost fringe cells
            ! right side is easy
            u0p1 = 0.0d0
            do j = 1,mshA%ngauss
              qvals = 0d0
              call shapefunction(2,mshA%xgauss(j),[-0.5d0,0.5d0],vec(n,:,eR),qvals,dqvals)
              u0p1 = u0p1 + SUM(qvals)*mshA%wgauss(j)
            enddo

            ! Get the left side from mesh B
            u0m1 = 0d0
            do j = 1,mshA%ngauss
              xg = mshA%xe(1,i)+ mshA%dx(i)*(mshA%xgauss(j)-0.5d0) ! x coord of gauss pt
              eloop1: do k = 1,mshB%nelem ! brute force loop through mesh B and find it's value
                xb1 = mshB%xe(1,k)
                xb2 = mshB%xe(2,k)
                if((xb1-xg)*(xb2-xg).le.0) then ! found cell containing gauss pt
                  dh = (xg-0.5d0*(xb1+xb2))/mshB%dx(k)
                  qvals = 0d0
                  call shapefunction(2,dh,[-0.5d0,0.5d0],mshB%q(n,:,k),qvals,dqvals)
                  u0m1 = u0m1 + sum(qvals)*mshA%wgauss(j)
                  exit eloop1 
                endif
              enddo eloop1 ! nelem 
            enddo ! ngauss

          else ! rightmost fringe cell
                  if(mshA%xe(1,i).gt.0d0) idebug = 1 
            ! left side is easy
            u0m1 = 0.0d0
            do j = 1,mshA%ngauss
              qvals = 0d0
              call shapefunction(mshA%nshp,mshA%xgauss(j),[-0.5d0,0.5d0],vec(n,:,eL),qvals,dqvals)
              u0m1 = u0m1 + SUM(qvals)*mshA%wgauss(j)
            enddo

            ! Get the right side from mesh B
            u0p1 = 0d0
            do j = 1,mshA%ngauss
              xg = mshA%xe(2,i)+mshA%dx(i)*(0.5d0 + mshA%xgauss(j)) ! x coord of gauss pt
              eloop2: do k = 1,mshB%nelem ! brute force loop through mesh B and find it's value
                xb1 = mshB%xe(1,k)
                xb2 = mshB%xe(2,k)
                if((xb1-xg)*(xb2-xg).le.0) then ! found cell containing gauss pt
                  dh = (xg-0.5d0*(xb1+xb2))/mshB%dx(k)
                  qvals = 0d0
                  call shapefunction(mshB%nshp,dh,[-0.5d0,0.5d0],mshB%q(n,:,k),qvals,dqvals)
                  u0p1 = u0p1 + sum(qvals)*mshA%wgauss(j)
                  exit eloop2
                endif
              enddo eloop2 ! nelem 
            enddo ! ngauss


          endif ! end of interior vs fringe cell case

          ! get the mean value of cell i
          u0   = 0.0d0
          do j = 1,mshA%ngauss
            qvals = 0d0
            call shapefunction(2,mshA%xgauss(j),[-0.5d0,0.5d0],vec(n,:,i),qvals,dqvals)
            u0   = u0   + SUM(qvals)*mshA%wgauss(j)
          enddo

          ! approximate un-limited end points 
          call shapefunction(mshA%nshp,-0.5d0,[-0.5d0,0.5d0],vec(n,:,i),qvals,dqvals)
          uL = sum(qvals);
          call shapefunction(mshA%nshp, 0.5d0,[-0.5d0,0.5d0],vec(n,:,i),qvals,dqvals)
          uR = sum(qvals); 

          ! get limited end point values
          call minmod2([u0-uL, u0-u0m1, u0p1 - u0],du,3,dx)
          uL2 = u0 - du 
          call minmod2([uR-u0, u0-u0m1, u0p1 - u0],du,3,dx)
          uR2 = u0 + du 

          ! if necessary, use MUSCL to limit q here
          diff = abs(uL-uL2) + abs(uR-uR2)
          tol = 1e-4
          if(diff.lt.tol) then ! keep orig solution
            vecout(n,:,i) = vec(n,:,i)
          else ! limit
            ! get linear slope of the element
            du0 = 0.0d0
            do j = 1,mshA%ngauss
              call shapefunction(mshA%nshp,mshA%xgauss(j),[-0.5d0,0.5d0],vec(n,:,i),qvals,dqvals)
              du0 = du0 + SUM(dqvals)*mshA%wgauss(j)/mshA%dx(i)
            enddo
  
            ! get the limited slope 
            call minmod2([du0,(u0p1-u0)/dx,(u0-u0m1)/dx],du,3,dx)

!          if(shptype.eq.'legendre') then 
              ! Now we have u (primary var) and need to reconstruct q (modal coeffs)
              ! int(Na Nb ub) = int(Na f(x,t))
              ! ub = M^-1 int(Na f(x,t))
              f = 0.0d0
              do j = 1,mshA%ngauss
                dh=mshA%xgauss(j)*dx
                tmp = u0 + dh*du ! limit to linear function
                qvals = 1d0
                call shapefunction(mshA%nshp,mshA%xgauss(j),[-0.5d0,0.5d0],qvals,qvals,dqvals)
                do k = 1,mshA%nshp
                  f(k) = f(k) + mshA%wgauss(j)*dx*tmp*qvals(k)
                enddo


              enddo
              ! solve the elementary system for u
              ! use the uncut mass here
              mass = 0.0
              do aa = 1,mshA%nshp
              do bb = 1,mshA%nshp
                index1 = (aa-1)*mshA%nshp+bb
                do k = 1,mshA%ngauss
                  mass(index1) = mass(index1) + mshA%shp(k,aa)*mshA%shp(k,bb)*mshA%wgauss(k)*mshA%dx(i)
                enddo
              enddo
              enddo

              call lu(mass,mshA%nshp,L,U)
              call forwprop(L,f,mshA%nshp,y)
              call backprop(U,y,mshA%nshp,vecout(n,:,i))
!              else ! Lagrange
             
!              endif ! elem type
          endif ! MUSCL limiter
        endif ! blanking 
      enddo ! loop over fields
    enddo ! loop over elements
  end subroutine genLimit2
!

  ! Pi^k Slope limiter from Cockburn (or Warburton pg 151) 
  subroutine genLimit(vec,vecout,msh)
    use code_types
    use pde
    use bases
    implicit none
    !
    type(mesh), intent(in) :: msh
    real*8,intent(in) :: vec(msh%nfields,msh%nshp,msh%nelem)
    real*8,intent(inout) :: vecout(msh%nfields,msh%nshp,msh%nelem)
    real*8 :: u0, u0m1, u0p1
    real*8 :: xg, xb1, xb2
    real*8 :: du0, du0m1, du0p1
    real*8 :: uL, uR, uL2, uR2, du
    real*8 :: qvals(msh%nshp), dqvals(msh%nshp)
    real*8 :: dh,tmp,f(msh%nshp),diff,tol,dx
    real*8,dimension(msh%nshp*msh%nshp) :: L,U
    real*8,dimension(msh%nshp) :: y
    integer :: n,i,j,k,eL,eR
    do i = 1,msh%nelem
      eL=msh%face(1,i) ! L neigh elem
      eR=msh%face(2,i) ! R neigh elem
      dx = msh%dx(i)
      do n = 1,msh%nfields
        if(maxval(msh%iblank(:,i)).eq.1) then ! ignore blanked cells
            ! Convert elements to linear and 
            ! get the cell averaged quantities u0 for elements i-1, i, and i+1
            u0m1 = 0.0d0
            u0p1 = 0.0d0
!            if(shptype.eq.'legendre') then 
              do j = 1,msh%ngauss
                ! need to get from other mesh if on boundary
                ! how does this work for overset?
                qvals = 0d0
                dqvals = 0d0
                call shapefunction(2,msh%xgauss(j),[-0.5d0,0.5d0],vec(n,:,eL),qvals,dqvals)
                u0m1 = u0m1 + SUM(qvals)*msh%wgauss(j) !jacobian and 1/area cancel out
                call shapefunction(2,msh%xgauss(j),[-0.5d0,0.5d0],vec(n,:,eR),qvals,dqvals)
                u0p1 = u0p1 + SUM(qvals)*msh%wgauss(j)
              enddo

          ! get the mean value of cell i
          u0   = 0.0d0
          do j = 1,msh%ngauss
            qvals = 0d0
            call shapefunction(2,msh%xgauss(j),[-0.5d0,0.5d0],vec(n,:,i),qvals,dqvals)
            u0   = u0   + SUM(qvals)*msh%wgauss(j)
          enddo

          ! approximate un-limited end points 
          call shapefunction(msh%nshp,-0.5d0,[-0.5d0,0.5d0],vec(n,:,i),qvals,dqvals)
          uL = sum(qvals);
          call shapefunction(msh%nshp, 0.5d0,[-0.5d0,0.5d0],vec(n,:,i),qvals,dqvals)
          uR = sum(qvals); 

          ! get limited end point values
          call minmod2([u0-uL, u0-u0m1, u0p1 - u0],du,3,dx)
          uL2 = u0 - du 
          call minmod2([uR-u0, u0-u0m1, u0p1 - u0],du,3,dx)
          uR2 = u0 + du 
  
          ! if necessary, use MUSCL to limit q here
          diff = abs(uL-uL2) + abs(uR-uR2)
          tol = 1e-4
          if(diff.lt.tol) then ! keep orig solution
            vecout(n,:,i) = vec(n,:,i)
          else ! limit
            ! get linear slope of the element
            du0 = 0.0d0
            do j = 1,msh%ngauss
              call shapefunction(msh%nshp,msh%xgauss(j),[-0.5d0,0.5d0],vec(n,:,i),qvals,dqvals)
              du0 = du0 + SUM(dqvals)*msh%wgauss(j)/msh%dx(i)
            enddo
  
            ! get the limited slope 
            call minmod2([du0,(u0p1-u0)/dx,(u0-u0m1)/dx],du,3,dx)

!          if(shptype.eq.'legendre') then 
              ! Now we have u (primary var) and need to reconstruct q (modal coeffs)
              ! int(Na Nb ub) = int(Na f(x,t))
              ! ub = M^-1 int(Na f(x,t))
              f = 0.0d0

              do j = 1,msh%ngauss
                dh=msh%xgauss(j)*dx
                tmp = u0 + dh*du ! limit to linear function
                qvals = 1d0
                call shapefunction(msh%nshp,msh%xgauss(j),[-0.5d0,0.5d0],qvals,qvals,dqvals)
                do k = 1,msh%nshp
                  f(k) = f(k) + msh%wgauss(j)*dx*tmp*qvals(k)
                enddo
              enddo

              ! solve the elementary system for u
              call lu(msh%mass(n,:,i),msh%nshp,L,U)
              call forwprop(L,f,msh%nshp,y)
              call backprop(U,y,msh%nshp,vecout(n,:,i))
!              else ! Lagrange
             
!              endif ! elem type
          endif ! MUSCL limiter
        endif ! blanking 
      enddo ! loop over fields
    enddo ! loop over elements
  end subroutine genLimit
!

end module slopelimiter
