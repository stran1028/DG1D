module code_types
  type mesh
     integer :: nfields
     integer :: nelem
     integer :: porder
     integer :: nnodes
     integer :: nfaces
     integer :: iperiodic
     integer :: nshp
     integer :: ngauss
     ! 
     integer, allocatable :: e2n(:,:)  !< element to end-node connectivity
     integer, allocatable :: face(:,:) !< element end-faces
     integer, allocatable :: iblank(:) !< iblank for end-nodes
     !
     real*8, allocatable :: dq(:,:,:)   !< q-variables (nfields,nshp,nelem)
     real*8, allocatable :: q(:,:,:)   !< intermediate q-variables (nfields,nshp,nelem)
     real*8, allocatable :: sol(:,:,:)   !< final q-variables (nfields,nshp,nelem)
     real*8, allocatable :: x(:)       !< end-node coords
     real*8, allocatable :: xe(:,:)    !< element end coords
     real*8, allocatable :: nres(:)    !< nodal resolution for overset
     real*8, allocatable :: rhs(:,:,:) !< rhs vector (nfield, nshp, nelem)
     real*8, allocatable :: rhsv(:,:,:) !< rhs vector (nfield, nshp, nelem)
     real*8, allocatable :: rhsf(:,:,:) !< rhs vector (nfield, nshp, nelem)
     real*8, allocatable :: mass(:,:,:,:) !< Mass matrix (nfield, nshp, nshp)
     real*8, allocatable :: dx(:)      !< element size
     real*8, allocatable :: wgauss(:)   !< weight of gauss pts (ngauss)
     real*8, allocatable :: xgauss(:)   !< location of gauss pts (ngauss)
     real*8, allocatable :: shp(:,:)     !< shape function values at qpts for parent element (ngauss,nshp)
     real*8, allocatable :: dshp(:,:)    !< shape function derivs for each shp fun at qpts (ngauss,nshp)
  end type mesh 
end module code_types
