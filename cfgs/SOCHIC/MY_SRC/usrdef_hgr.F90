MODULE usrdef_hgr
   !!======================================================================
   !!                     ***  MODULE usrdef_hgr   ***
   !!
   !!                     ===  GYRE configuration  ===
   !!
   !! User defined :   mesh and Coriolis parameter of a user configuration
   !!======================================================================
   !! History :  4.0 ! 2016-03  (S. Flavoni) 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_hgr   : initialize the horizontal mesh 
   !!----------------------------------------------------------------------
   USE dom_oce  , ONLY: nimpp, njmpp       ! ocean space and time domain
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE usrdef_nam     !
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_hgr   ! called in domhgr.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_hgr.F90 10069 2018-08-28 14:12:24Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_hgr( plamt , plamu , plamv  , plamf  ,   &   ! geographic position (required)
      &                    pphit , pphiu , pphiv  , pphif  ,   &   !
      &                    kff   , pff_f , pff_t  ,            &   ! Coriolis parameter  (if domain not on the sphere)
      &                    pe1t  , pe1u  , pe1v   , pe1f   ,   &   ! scale factors       (required)
      &                    pe2t  , pe2u  , pe2v   , pe2f   ,   &   !
      &                    ke1e2u_v      , pe1e2u , pe1e2v     )   ! u- & v-surfaces (if gridsize reduction is used in strait(s))
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE usr_def_hgr  ***
      !!
      !! ** Purpose :   user defined mesh and Coriolis parameter
      !!
      !! ** Method  :   set all intent(out) argument to a proper value
      !!
      !!                Here GYRE configuration :
      !!          Rectangular mid-latitude domain 
      !!          - with axes rotated by 45 degrees
      !!          - a constant horizontal resolution of 106 km 
      !!          - on a beta-plane
      !!
      !! ** Action  : - define longitude & latitude of t-, u-, v- and f-points (in degrees) 
      !!              - define coriolis parameter at f-point if the domain in not on the sphere (on beta-plane)
      !!              - define i- & j-scale factors at t-, u-, v- and f-points (in meters)
      !!              - define u- & v-surfaces (if gridsize reduction is used in some straits) (in m2)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   plamt, plamu, plamv, plamf   ! longitude outputs                     [degrees]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pphit, pphiu, pphiv, pphif   ! latitude outputs                      [degrees]
      INTEGER                 , INTENT(out) ::   kff                          ! =1 Coriolis parameter computed here, =0 otherwise
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pff_f, pff_t                 ! Coriolis factor at f-point                [1/s]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1t, pe1u, pe1v, pe1f       ! i-scale factors                             [m]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe2t, pe2u, pe2v, pe2f       ! j-scale factors                             [m]
      INTEGER                 , INTENT(out) ::   ke1e2u_v                     ! =1 u- & v-surfaces computed here, =0 otherwise 
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1e2u, pe1e2v               ! u- & v-surfaces (if reduction in strait)   [m2]
      !
      INTEGER  ::   ji, jj               ! dummy loop indices
      REAL(wp) ::   x0, zcos_alpha, zim1 , zjm1 , ze1  , ze1deg, zf0 ! local scalars
      REAL(wp) ::   zphi, y0, zsin_alpha, zim05, zjm05, zbeta, znorme      !   -      -
      !!-------------------------------------------------------------------------------
      !
      !     !==  beta-plane with regular grid-spacing and rotated domain ==!  (GYRE configuration)
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_hgr : GYRE configuration (beta-plane with rotated regular grid-spacing)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) 'grid hgr dy: ', dy
      !
      !
      !                       !==  grid point position  ==!
      !
      zphi = -60._wp
      !
      !dx = 2000._wp   ! gridspacing in meters
      !dy = 2000._wp   ! gridspacing in meters
      !
      x0 = -(jpiglo)/2 * dx * 1.e-3 ! x0
      y0 = -(jpjglo)/2 * dy * 1.e-3 ! y0

!#if defined key_agrif
!      ! ! Upper left longitude and latitude from parent:
!      IF (.NOT.Agrif_root()) THEN
!         zlam = zlam + Agrif_irhox() * REAL(Agrif_Parent(jpjglo)-2 , wp) * ze1deg * zcos_alpha  &
!                   &   + ( Agrif_Ix()*Agrif_irhox()-(0.5_wp+nbghostcells)) * ze1deg * zcos_alpha  &
!                   &   + ( Agrif_Iy()*Agrif_irhoy()-(0.5_wp+nbghostcells)) * ze1deg * zsin_alpha
!         zphi = zphi + Agrif_irhoy() * REAL(Agrif_Parent(jpjglo)-2 , wp) * ze1deg * zsin_alpha  &
!                   &   - ( Agrif_Ix()*Agrif_irhox()-nbghostcells )         * ze1deg * zsin_alpha  &
!                   &   + ( Agrif_Iy()*Agrif_irhoy()-nbghostcells )         * ze1deg * zcos_alpha
!      ENDIF 
!#endif
!      !   
      IF( nprint==1 .AND. lwp )   THEN
         WRITE(numout,*) 'dx', dx, 'x0', x0, 'y0', y0
      ENDIF
      !   
      DO jj = 1, jpj 
         DO ji = 1, jpi 
            zim1 = FLOAT( ji + nimpp - 1 ) - 1.   ; zim05 = FLOAT( ji + nimpp - 1 ) - 0.5_wp 
            zjm1 = FLOAT( jj + njmpp - 1 ) - 1.   ; zjm05 = FLOAT( jj + njmpp - 1 ) - 0.5_wp 
            !   
            !glamt(i,j) longitude at T-point
            !gphit(i,j) latitude at T-point  
            plamt(ji,jj) = x0 + zim1 * dx * 1.e-3
            pphit(ji,jj) = y0 + zjm1 * dy * 1.e-3
            !   
            !glamu(i,j) longitude at U-point
            !gphiu(i,j) latitude at U-point
            plamu(ji,jj) = x0 + zim05 *  dx * 1.e-3
            pphiu(ji,jj) = pphit(ji,jj)
            !   
            !glamv(i,j) longitude at V-point
            !gphiv(i,j) latitude at V-point
            plamv(ji,jj) = plamt(ji,jj)
            pphiv(ji,jj) = y0 + zjm05 * dy * 1.e-3
            !
            !glamf(i,j) longitude at F-point
            !gphif(i,j) latitude at F-point 
            plamf(ji,jj) = plamu(ji,jj)
            pphif(ji,jj) = pphiv(ji,jj)
         END DO
      END DO
      !
      !                       !== Horizontal scale factors ==! (in meters)
      !                     
      !                                         ! constant grid spacing
      pe1t(:,:) =  dx     ;      pe2t(:,:) = dy
      pe1u(:,:) =  dx     ;      pe2u(:,:) = dy
      pe1v(:,:) =  dx     ;      pe2v(:,:) = dy
      pe1f(:,:) =  dx     ;      pe2f(:,:) = dy
      !
      !                                         ! NO reduction of grid size in some straits 
      ke1e2u_v = 0        !    ==>> u_ & v_surfaces will be computed in dom_ghr routine
      pe1e2u(:,:) = 0._wp                       !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v(:,:) = 0._wp                       !             require an initialization of INTENT(out) arguments
      !
      !
      !                       !==  Coriolis parameter  ==!
      kff = 1                                            !  indicate not to compute ff afterward
      !
      zbeta = 2._wp * omega * COS( rad * zphi ) / ra       ! beta at latitude zphi1
      zf0   = 2._wp * omega * SIN( rad * zphi )            !  compute f0 1st point south
      !
      pff_f(:,:) =  zf0 + zbeta * pphif(:,:) * 1.e+3 ! f = f0 +beta* y ( y=0 at south)
      pff_t(:,:) =  zf0 + zbeta * pphit(:,:) * 1.e+3 ! f = f0 +beta* y ( y=0 at south)
      !
      IF(lwp) WRITE(numout,*) '                           beta-plane used. beta = ', zbeta, ' 1/(s.m)'
      IF(lwp) WRITE(numout,*) '                           f0 = ', pff_t(0,:), ' 1/(s.m)'
      !
   END SUBROUTINE usr_def_hgr

   !!======================================================================
END MODULE usrdef_hgr
