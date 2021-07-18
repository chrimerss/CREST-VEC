MODULE basinUH_module

USE nrtype
USE public_var
USE dataTypes,          only : STRFLX         ! fluxes in each reach

implicit none

private

public::IRF_route_basin, IRF_route_basin_subsurf

CONTAINS

 ! ---------------------------------------------------------------------------------------
 ! Public subroutine main driver for basin routing
 ! ---------------------------------------------------------------------------------------
 SUBROUTINE IRF_route_basin(iens,          & ! input: ensemble index
                            reachRunoff,   & ! input: instantaneous reach runoff volume from hru (m3/s)
                            flag_subsurf,  & ! input: flag whether activate subsurface route
                            RCHFLX_out,    & ! inout: reach flux data structure
                            ierr, message, & ! output: error control
                            ixSubRch)        ! optional input: subset of reach indices to be processed
 implicit none
 ! input
 integer(i4b), intent(in)                 :: iens            ! ith ensemble
 real(dp)    , intent(in)                 :: reachRunoff(:)  ! instantaneous reach runoff volumen (m3/s)
 LOGICAL, intent(in)                      :: flag_subsurf    ! whether activates subsurface route
 ! inout
 type(STRFLX), intent(inout), allocatable :: RCHFLX_out(:,:) ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! output
 integer(I4B), intent(out)                :: ierr            ! error code
 character(*), intent(out)                :: message         ! error message
 ! input (optional)
 integer(i4b), intent(in),   optional     :: ixSubRch(:)     ! subset of reach indices to be processed
 ! local variables
 integer(i4b)                             :: nSeg            ! number of reaches to be processed
 integer(i4b)                             :: iSeg            ! reach loop indix
 logical(lgt),               allocatable  :: doRoute(:)      ! logical to indicate which reaches are processed
 character(len=strLen)                    :: cmessage        ! error message from subroutines

 ierr=0; message='IRF_route_basin/'

 nSeg = size(RCHFLX_out(iens,:))

 allocate(doRoute(nSeg), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'unable to allocate space for [doRoute]'; return; endif

 ! if a subset of reaches is processed
 if (present(ixSubRch))then
  doRoute(:) = .false.
  doRoute(ixSubRch) = .true. ! only subset of reaches are on
 else
  doRoute(:) = .true.
 endif

!$OMP PARALLEL DO schedule(dynamic,1)   &
!$OMP          private(iSeg)            & ! loop index
!$OMP          private(ierr, cmessage)  & ! private for a given thread
!$OMP          shared(doRoute)          & ! data array shared
!$OMP          shared(RCHFLX_out)       & ! data structure shared
!$OMP          shared(iEns)             & ! indices shared
!$OMP          firstprivate(nSeg)
 do iSeg=1,nSeg

  if (.not. doRoute(iSeg)) cycle

  RCHFLX_out(iens,iSeg)%BASIN_QI = reachRunoff(iSeg)
  call hru_irf(iEns, iSeg, flag_subsurf, RCHFLX_out, ierr, cmessage)
!  f(ierr/=0)then; ixmessage(iSeg)=trim(message)//trim(cmessage); exit; endif

 end do
!$OMP END PARALLEL DO

 END SUBROUTINE IRF_route_basin


 ! ---------------------------------------------------------------------------------------
 ! Public subroutine main driver for basin subsurface routing
 ! ---------------------------------------------------------------------------------------
 SUBROUTINE IRF_route_basin_subsurf(iens,          & ! input: ensemble index
                                    reachSubRunoff,   & ! input: instantaneous reach runoff volume from hru (m3/s)
                                    flag_subsurf,  & ! input: flag whether activate subsurface route                                    
                                    RCHFLX_SUB,    & ! inout: reach flux data structure
                                    ierr, message, & ! output: error control
                                    ixSubRch)        ! optional input: subset of reach indices to be processed
 implicit none
 ! input
 integer(i4b), intent(in)                 :: iens            ! ith ensemble
 real(dp)    , intent(in)                 :: reachSubRunoff(:)  ! instantaneous reach runoff volumen (m3/s)
 LOGICAL, intent(in)                      :: flag_subsurf    ! whether activates subsurface route
 ! inout
 type(STRFLX), intent(inout), allocatable :: RCHFLX_SUB(:,:) ! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! output
 integer(I4B), intent(out)                :: ierr            ! error code
 character(*), intent(out)                :: message         ! error message
 ! input (optional)
 integer(i4b), intent(in),   optional     :: ixSubRch(:)     ! subset of reach indices to be processed
 ! local variables
 integer(i4b)                             :: nSeg            ! number of reaches to be processed
 integer(i4b)                             :: iSeg            ! reach loop index
 logical(lgt),               allocatable  :: doRoute(:)      ! logical to indicate which reaches are processed
 character(len=strLen)                    :: cmessage        ! error message from subroutines

 ierr=0; message='IRF_route_basin_subsurf/'

 nSeg = size(RCHFLX_SUB(iens,:))

 allocate(doRoute(nSeg), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'unable to allocate space for [doRoute]'; return; endif


 ! if a subset of reaches is processed
 if (present(ixSubRch))then
  doRoute(:) = .false.
  doRoute(ixSubRch) = .true. ! only subset of reaches are on
 else
  doRoute(:) = .true.
 endif

!$OMP PARALLEL DO schedule(dynamic,1)   &
!$OMP          private(iSeg)            & ! loop index
!$OMP          private(ierr, cmessage)  & ! private for a given thread
!$OMP          shared(doRoute)          & ! data array shared
!$OMP          shared(RCHFLX_out)       & ! data structure shared
!$OMP          shared(iEns)             & ! indices shared
!$OMP          firstprivate(nSeg)
 do iSeg=1,nSeg

  if (.not. doRoute(iSeg)) cycle

  RCHFLX_SUB(iens,iSeg)%BASIN_QI = reachSubRunoff(iSeg)
  call hru_irf(iEns, iSeg, flag_subsurf, RCHFLX_SUB, ierr, cmessage)
  ! if(ierr/=0)then; ixmessage(iSeg)=trim(message)//trim(cmessage); exit; endif
  ! RCHFLX_out(iens,iSeg)%TAKE= RCHSUBSURF(iens,iSeg)%BASIN_QR(0)
 end do
 !$OMP END PARALLEL DO

 END SUBROUTINE IRF_route_basin_subsurf

 ! *********************************************************************
 ! subroutine: perform one basin UH routing to a iSeg reach at one time
 ! *********************************************************************
 subroutine hru_irf(iens,         &    ! input: index of runoff ensemble to be processed
                    iSeg,         &    ! input: index of runoff ensemble to be processed
                    flag_subsurf, &    ! input: whether it is in subsurface routing
                    RCHFLX_out,   &    ! inout: reach flux data structure
                    ierr, message)     ! output: error control
 ! External modules
 USE globalData,        ONLY : FRAC_FUTURE, FRAC_FUTURE_SUB
 implicit none
 ! Input
 INTEGER(I4B), intent(IN)                 :: iEns           ! runoff ensemble to be routed
 INTEGER(I4B), intent(IN)                 :: iSeg           ! segment where routing is performed
 LOGICAL                                  :: flag_subsurf   ! whether it is in subsurface routing
 ! inout
 TYPE(STRFLX), intent(inout), allocatable :: RCHFLX_out(:,:)! Reach fluxes (ensembles, space [reaches]) for decomposed domains
 ! Output
 integer(i4b), intent(out)                :: ierr           ! error code
 character(*), intent(out)                :: message        ! error message
 ! Local variables to
 INTEGER(I4B)                             :: ntdh           ! number of time steps in IRF
 character(len=strLen)                    :: cmessage       ! error message from subroutine

 ! initialize error control
 ierr=0; message='hru_irf/'

 ! initialize the first time step q future
  if (.not.allocated(RCHFLX_out(iens,iSeg)%QFUTURE))then
   if (flag_subsurf) then
    ntdh = size(FRAC_FUTURE_SUB)
   else
    ntdh= size(FRAC_FUTURE)
   endif
   allocate(RCHFLX_out(iens,iSeg)%QFUTURE(ntdh), stat=ierr)
   if(ierr/=0)then; message=trim(message)//'unable to allocate space for RCHFLX_out(iens,segIndex)%QFUTURE'; return; endif


   RCHFLX_out(iens,iSeg)%QFUTURE(:) = 0._dp

  end if

  ! perform river network UH routing
  if (flag_subsurf) then
    call irf_conv(FRAC_FUTURE_SUB,                     &    ! input: unit hydrograph
                  RCHFLX_out(iens,iSeg)%BASIN_QI,  &    ! input: upstream fluxes
                  RCHFLX_out(iens,iSeg)%QFUTURE,   &    ! inout: updated q future time series
                  RCHFLX_out(iens,iSeg)%BASIN_QR,  &    ! inout: updated fluxes at reach
                ierr, message)                            ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  else
    call irf_conv(FRAC_FUTURE,                     &    ! input: unit hydrograph
                  RCHFLX_out(iens,iSeg)%BASIN_QI,  &    ! input: upstream fluxes
                  RCHFLX_out(iens,iSeg)%QFUTURE,   &    ! inout: updated q future time series
                  RCHFLX_out(iens,iSeg)%BASIN_QR,  &    ! inout: updated fluxes at reach
                ierr, message)                            ! output: error control
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

 end subroutine hru_irf


 ! ---------------------------------------------------------------------------------------
 ! Private subroutine: Perform UH convolutions
 ! ---------------------------------------------------------------------------------------
 SUBROUTINE irf_conv(uh,        &    ! input: normalized unit hydrograph
                     inq,       &    ! input: instantaneous runoff
                     qfuture,   &    ! inout: convoluted runoff including future time step
                     delayq,    &    ! inout: delayed runoff to segment at a current and previous time step
                     ierr, message)
  implicit none
  ! input
  real(dp),             intent(in)     :: uh(:)         ! normalized unit hydrograph
  real(dp),             intent(in)     :: inq           ! basin instantaneous runoff
  real(dp),             intent(inout)  :: qfuture(:)    ! convoluted runoff including future time steps
  real(dp),             intent(inout)  :: delayq(0:1)   ! delayed runoff to a segment at a current and previous time step
  ! output
  integer(I4B),         intent(out)    :: ierr          ! error code
  character(*),         intent(out)    :: message       ! error message
  ! local variables
  integer(i4b)                         :: itdh          ! index loop for basin, time, respectively
  integer(i4b)                         :: ntdh          ! number of time step for future flow

  ! initialize error control
  ierr=0; message='irf_conv/'

  ! route streamflow through the basin
  ntdh = size(qfuture)

  ! place a fraction of runoff in future time steps and add to current state of q in the future
  do itdh=1,ntdh
   qfuture(itdh) = qfuture(itdh) + uh(itdh)*inq
  end do

  ! save the routed runoff
  delayq(0) = delayq(1)  ! (save the runoff from the previous time step)
  delayq(1) = qfuture(1)

  ! move array back
  do itdh=2,ntdh
   qfuture(itdh-1) = qfuture(itdh)
  end do
  qfuture(ntdh)    = 0._dp

  END SUBROUTINE irf_conv

END MODULE basinUH_module
