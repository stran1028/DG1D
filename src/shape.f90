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
        endif
      elseif(shptype.eq.'legendre') then 
        ! Evaluate legendre polynomial
!        write(*,*) 'LEGENDRE POLY'
        do i = 1,nshp
          do j = 1,i
            qval(i) = qval(i) + q(i)*legcoef(i,j)*zi**(j-1d0)
            dqval(i) = dqval(i) + q(i)*legderv(i,j)*zi**(j-1d0)
!            if(i.eq.nshp) then 
!             write(*,*) 'j = ',j
!             write(*,*) '  coeff1 = ',legcoef(i,j)
!             write(*,*) '  coeff2 = ',legderv(i,j)
!             write(*,*) '  power = ',j-1d0
!             write(*,*) '  q = ',q(i)
!             write(*,*) '  zi = ',zi
!             write(*,*) '  tmp = ',(2d0*zi)
!             write(*,*) '  tmp1 = ',(2d0*zi)**(j-1d0)
!             write(*,*) '  tmp2 = ',q(i)*legcoef(i,j)
!             write(*,*) '  qvaltmp = ',qval(i)
!            endif
          enddo
!          write(*,*) 'i, qval = ',i,qval(i)
        enddo
        
!        test(1) = 1d0*q(1)
!        dtest(1) = 0d0*q(1)
!        test(2) = 2d0*zi*q(2) 
!        dtest(2) = 2d0*q(2)
!        do i = 1,2
!        if((abs(test(i)-qval(i)).gt.1e-8).or.(abs(dtest(i)-dqval(i)).gt.1e-8)) then 
!          write(*,*) 'error: ',test(i),qval(i)
!          write(*,*) 'derror: ',dtest(i),dqval(i)
!          call exit(1)
!        endif
!        enddo
      !  qval(3) = 0.5d0*(3d0*zi*zi-1d0) *q(3)
!        dqval(3) = 3d0*zi*q(3)
      endif
      qout = qval(1:nshp)
      dqout = dqval(1:nshp)

    end subroutine shapefunction
end module bases
