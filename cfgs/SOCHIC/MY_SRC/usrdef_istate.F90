MODULE usrdef_istate
   !!======================================================================
   !!                   ***  MODULE  usrdef_istate   ***
   !!
   !!                     ===  GYRE configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  4.0 ! 2016-03  (S. Flavoni) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   !USE dom_oce  , ONLY: nimpp, njmpp       ! ocean space and time domain
   USE usrdef_nam     !
!   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
!   USE sbc_oce        ! Surface boundary condition: ocean fields

   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate   ! called in istate.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_istate.F90 10069 2018-08-28 14:12:24Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here GYRE configuration example : (double gyre with rotated domain)
      !!
      !! ** Method  : - set temprature field
      !!              - set salinity   field
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height
      !
      INTEGER  :: ji, jj, jk  ! dummy loop indices
      !INTEGER, DIMENSION(ji,jj,jk) :: seed
      REAL(wp) :: rand       ! random number
      REAL(wp) :: maxlat      ! grid spacing in y


      ! ---------------------------- !
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : analytical definition of initial state '
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   Ocean at rest, with an horizontally uniform T and S profiles'
      IF(lwp) WRITE(numout,*) 'grid istate dy: ', dy
      !
!      pu  (:,:,:) = 0._wp
!      pv  (:,:,:) = 0._wp
      pssh(:,:) = 0._wp
      !

!      DO jj = 1, jpj
!         DO ji = 1, jpi
!            CALL RANDOM_NUMBER(rand)
!            pssh(ji,jj)  = ( rand - 0.5_wp ) * 1.e-1
!         END DO
!      END DO

      CALL RANDOM_SEED()
      maxlat = jpjglo/2 * dy * 1.e-3 ! x0
      DO jk = 1, jpk             ! horizontally uniform T & S profiles
       DO jj = 1, jpj
        DO ji = 1, jpi
         CALL RANDOM_Number(rand)
         pu(ji,jj,jk) = 0.2 + ( rand - 0.5_wp ) * 1.e-3
         pv(ji,jj,jk) = ( rand - 0.5_wp ) * 1.e-3
         pts(ji,jj,jk,jp_tem) =                                         &
    &     ( 1.5 + (0.5  - 1.0 * (1 + gphit(ji,jj) / (2 * maxlat) ))      &
    &     * 0.5 *  ( TANH ((pdept(ji,jj,jk) - 100) / 50 ) - 1 )         &
    &                                            ) * ptmask(ji,jj,jk)

         pts(ji,jj,jk,jp_sal) =                                         &
    &     ( 34.5 + (0.2 + 2.0 * (1 + gphit(ji,jj) / (2 * maxlat) ))     &
    &     *  0.5 *  ( TANH ((pdept(ji,jj,jk) -100) / 50 ) - 1 )         &
    &                  ) * ptmask(ji,jj,jk)
!         pts(ji,jj,jk,jp_tem) = 5._wp
!         pts(ji,jj,jk,jp_sal) = 34.5_wp
        END DO
       END DO
      END DO
      !   
   END SUBROUTINE usr_def_istate

   !!======================================================================
END MODULE usrdef_istate
