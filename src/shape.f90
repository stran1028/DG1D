module bases
  character*20 :: shptype

contains
    subroutine setshp(val)
      character*(*) :: val
      shptype = trim(adjustl(val))
      write(*,*) 'shp = ',shptype
      if((shptype.ne.'lagrange').and.(shptype.ne.'legendre')) then
        write(*,*) 'That shape function is not supported. Exiting.'
        call exit(1)
      endif
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
      real*8 :: dx,xloc,xc,zi,qval(5),dqval(5)
      !
      ! Transform element to go from [-.5 .5]
      dx = xlim(2)-xlim(1)
      xc = 0.5d0*(xlim(2)+xlim(1))
      zi = (x-xc)/dx !<+ wrong
      !write(*,*) 'x,dx,xc,zi',x,dx,xc,zi
      !
      qval = 0.0d0
      if(shptype.eq.'lagrange') then 
        qval(1)  = q(1)*(0.5d0-zi)
        dqval(1) = -1d0*q(1)
        qval(2)  = q(2)*(0.5d0+zi)
        dqval(2) = 1d0*q(2)
        qout = qval(1:nshp)
        dqout = dqval(1:nshp)
      elseif(shptype.eq.'legendre') then 
        qval(1) = 1d0*q(1)
        dqval(1) = 0d0*q(1)
        qval(2) = zi*q(2) 
        dqval(2) = 1d0*q(2)
       ! write(*,*) 'qval = ',qval
       ! write(*,*) 'dqval = ',dqval
!        qval(3) = 0.5d0*(3d0*zi*zi-1d0) *q(3)
!        dqval(3) = 3d0*zi*q(3)
        qout = qval(1:nshp)
        dqout = dqval(1:nshp)
      endif

    end subroutine shapefunction
end module bases
