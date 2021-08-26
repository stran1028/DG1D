module bases
  character*20 :: shptype
  real*8 :: legcoef(11,11)
  real*8 :: legderv(11,11)

contains
    subroutine setshp(val)
      character*(*) :: val
      integer :: i,j
      shptype = trim(adjustl(val))
      if((shptype.ne.'lagrange').and.(shptype.ne.'legendre')) then
        write(*,*) 'That shape function is not supported. Exiting.'
        call exit(1)
      endif

      ! Load in legendre coefficients from Beatrice's files
      open(1,file = 'LegendreCoeffs1.dat',status='old')
      do i = 1,11
        read(1,*) legcoef(i,:)
      enddo
      close(1)

      open(1,file = 'LegendreCoeffs2.dat',status='old')
      do i = 1,11
        read(1,*) legderv(i,:)
      enddo
      close(1)

    end subroutine setshp
    !
    subroutine shapefunction(nshp,x,xlim,q,qout,dqout)
      implicit none
      !
      integer, intent(in) :: nshp
      real*8, intent(in) :: xlim(2),x,q(nshp)
      real*8, intent(inout) :: qout(nshp),dqout(nshp)
      !
      integer:: i,j,k
      real*8 :: dx,xloc,xc,zi,qval(11),dqval(11)
      real*8 :: test(2),dtest(2)
      !
      ! Transform element to go from [-.5 .5]
      dx = xlim(2)-xlim(1)
      xc = 0.5d0*(xlim(2)+xlim(1))
      zi = (x-xc)/dx 
      !
      qval = 0d0
      dqval = 0d0
      if(shptype.eq.'lagrange') then 
        if(nshp.eq.1) then 
          qval(1) = 1d0*q(1)
          dqval(1) = 0d0*q(1)
        elseif(nshp.eq.2) then 
          qval(1)  = q(1)*(0.5d0-zi)
          dqval(1) = -1d0*q(1)
          qval(2)  = q(2)*(0.5d0+zi)
          dqval(2) = 1d0*q(2)
        elseif(nshp.eq.3) then 
          qval(1)  = q(1)*(2d0*zi**2d0 - zi)
          dqval(1) = q(1)*(4d0*zi-1d0)
          qval(2)  = q(2)*(-4d0*zi**2d0 + 1d0)
          dqval(2) = q(2)*(-8d0*zi)
          qval(3)  = q(3)*(2d0*zi**2d0 + zi)
          dqval(3) = q(3)*(4d0*zi+1d0) 
        elseif(nshp.eq.4) then 
          qval(1)  = q(1)*(-4.5d0*zi**3d0 + 2.25*zi**2d0 + 0.125*zi-0.0625)
          dqval(1) = q(1)*(-13.5d0*zi**2d0 + 4.5d0*zi + 0.125)
          qval(2)  = q(2)*(13.5d0*zi**3d0 - 2.25d0*zi**2d0 - 3.375*zi + 0.5625)
          dqval(2) = q(2)*( 40.5d0*zi**2d0 - 4.5d0*zi - 3.375)
          qval(3)  = q(3)*(-13.5d0*zi**3d0 - 2.25d0*zi**2d0 + 3.375*zi + 0.5625)
          dqval(3) = q(3)*(-40.5d0*zi**2d0 - 4.5d0*zi + 3.375)
          qval(4)  = q(4)*( 4.5d0*zi**3d0 + 2.25*zi**2d0 - 0.125*zi-0.0625)
          dqval(4) = q(4)*( 13.5d0*zi**2d0 + 4.5d0*zi - 0.125)
        elseif(nshp.eq.5) then 
          qval(1) = q(1)*(32d0/3d0*zi**4d0 - 16d0/3d0*zi**3d0 -2d0/3d0*zi**2d0 &
                        + 1d0/3d0*zi)
          qval(2) = q(2)*(-128d0/3d0*zi**4d0 + 32d0/3d0*zi**3d0 &
                          + 32d0/3d0*zi**2d0 - 8d0/3d0*zi)
          qval(3) = q(3)*(64d0*zi**4d0 - 20d0*zi**2d0 + 1d0)
          qval(4) = q(4)*(-128d0/3d0*zi**4d0 - 32d0/3d0*zi**3d0 &
                          + 32d0/3d0*zi**2d0 + 8d0/3d0*zi)
          qval(5) = q(5)*(32d0/3d0*zi**4d0 + 16d0/3d0*zi**3d0 -2d0/3d0*zi**2d0 &
                        - 1d0/3d0*zi)
          dqval(1) = q(1)*(128d0/3d0*zi**3d0 - 16d0*zi**2d0-4d0/3d0*zi+1d0/3d0 )
          dqval(2) = q(2)*(-512d0/3d0*zi**3d0+32d0*zi**2d0 + 64d0/3d0*zi-8d0/3d0)
          dqval(3) = q(3)*(256d0*zi**3d0-40d0*zi) 
          dqval(4) = q(4)*(-512d0/3d0*zi**3d0-32d0*zi**2d0 + 64d0/3d0*zi+8d0/3d0)
          dqval(5) = q(5)*(128d0/3d0*zi**3d0 + 16d0*zi**2d0-4d0/3d0*zi-1d0/3d0 )
        elseif(nshp.eq.6) then 
          qval(1) = q(1)*( -625d0/24d0*zi**5d0 + 625d0/48d0*zi**4d0 &
                           +125d0/48d0*zi**3d0 - 125d0/96d0*zi**2d0 &
                           - 3d0/128d0*zi      + 3d0/256d0)
          qval(2) = q(2)*( 3125d0/24d0*zi**5d0 - 625d0/16d0*zi**4d0 & 
                          -1625d0/48d0*zi**3d0 + 325d0/32d0*zi**2d0 &
                          +125d0/384d0*zi      - 25d0/256d0)
          qval(3) = q(3)*(-3125d0/12d0*zi**5d0 + 625d0/24d0*zi**4d0 &
                          +2125d0/24d0*zi**3d0 - 425d0/48d0*zi**2d0 &
                          -375d0/64d0*zi       + 75d0/128d0)
          qval(4) = q(4)*( 3125d0/12d0*zi**5d0 + 625d0/24d0*zi**4d0 &
                          -2125d0/24d0*zi**3d0 - 425d0/48d0*zi**2d0 &
                          +375d0/64d0*zi       + 75d0/128d0)
          qval(5) = q(5)*(-3125d0/24d0*zi**5d0 - 625d0/16d0*zi**4d0 & 
                          +1625d0/48d0*zi**3d0 + 325d0/32d0*zi**2d0 &
                          -125d0/384d0*zi      - 25d0/256d0)
          qval(6) = q(6)*(  625d0/24d0*zi**5d0 + 625d0/48d0*zi**4d0 &
                           -125d0/48d0*zi**3d0 - 125d0/96d0*zi**2d0 &
                           + 3d0/128d0*zi      + 3d0/256d0)

          dqval(1) = q(1)*(-3125d0/24d0*zi**4d0+625d0/12d0*zi**3d0 &
                            +125d0/16d0*zi**2d0-125d0/48d0*zi&
                            -3d0/128d0)
          dqval(2) = q(2)*(15625d0/24d0*zi**4d0-625d0/4d0*zi**3d0 &
                           -1625d0/16d0*zi**2d0+325d0/16d0*zi&
                           +125d0/384d0)
          dqval(3) = q(3)*(-15625d0/12d0*zi**4d0+625d0/6d0*zi**3d0 &
                           +2125d0/8d0*zi**2d0 -425d0/24d0*zi &
                           - 375d0/64d0)
          dqval(4) = q(4)*(15625d0/12d0*zi**4d0+625d0/6d0*zi**3d0 &
                          -2125d0/8d0*zi**2d0 -425d0/24d0*zi &
                           + 375d0/64d0)
          dqval(5) = q(5)*(-15625d0/24d0*zi**4d0-625d0/4d0*zi**3d0 &
                           +1625d0/16d0*zi**2d0+325d0/16d0*zi &
                           -125d0/384d0)
          dqval(6) = q(6)*(3125d0/24d0*zi**4d0+625d0/12d0*zi**3d0 &
                           -125d0/16d0*zi**2d0 - 125d0/48d0*zi &
                           +3d0/128d0)
        elseif(nshp.gt.6) then 
          write(*,*) 'Shape functions not yet implemented. Exiting'
          call exit(1)
        endif
      elseif(shptype.eq.'legendre') then 
        ! Evaluate legendre polynomial
        if(nshp.gt.10) then 
          write(*,*) 'Shape functions not yet implemented. Exiting'
          call exit(1)
        endif
        do i = 1,nshp
          do j = 1,i
            qval(i) = qval(i) + q(i)*legcoef(i,j)*zi**(j-1d0)
            dqval(i) = dqval(i) + q(i)*legderv(i,j)*zi**(j-1d0)
          enddo
        enddo
      endif
      qout = qval(1:nshp)
      dqout = dqval(1:nshp)

    end subroutine shapefunction
end module bases
