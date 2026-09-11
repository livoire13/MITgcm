#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: EXTERNAL_FIELDS_LOAD
C     !INTERFACE:
      SUBROUTINE EXTERNAL_FIELDS_LOAD( myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE EXTERNAL_FIELDS_LOAD
C     | o Control reading of fields from external source.
C     *==========================================================*
C     | External source field loading routine.
C     | This routine is called every time we want to
C     | load a a set of external fields. The routine decides
C     | which fields to load and then reads them in.
C     | This routine needs to be customised for particular
C     | experiments.
C     | Notes
C     | =====
C     | Two-dimensional and three-dimensional I/O are handled in
C     | the following way under MITgcmUV. A master thread
C     | performs I/O using system calls. This threads reads data
C     | into a temporary buffer. At present the buffer is loaded
C     | with the entire model domain. This is probably OK for now
C     | Each thread then copies data from the buffer to the
C     | region of the proper array it is responsible for.
C     | =====
C     | Conversion of flux fields are described in FFIELDS.h
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "FFIELDS.h"
#include "GRID.h"
#include "DYNVARS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime :: Simulation time
C     myIter :: Simulation timestep number
C     myThid :: Thread no. that called this routine.
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

#ifndef EXCLUDE_FFIELDS_LOAD

C     !LOCAL VARIABLES:
C     === Local arrays ===
C     aWght, bWght :: Interpolation weights
      INTEGER bi, bj, i, j, k
      INTEGER intimeP, intime0, intime1
      _RL aWght, bWght
CHL   !LOCAL VARIABLES FOR SGS FORCING:
      INTEGER sgsIntimeP, sgsIntime0, sgsIntime1
      _RL sgsAWght, sgsBWght
      LOGICAL sgsExternalForcing
CEOP

CHL   Set sgsExternalForcing boolean
      IF ( sgsThetaFile .NE. ' ' .OR.
     &     sgsSaltFile  .NE. ' ' ) THEN
        sgsExternalForcing = .TRUE.
      ENDIF

      IF ( periodicExternalForcing ) THEN

C--   First call requires that we initialize everything to zero for safety
cph   has been shifted to ini_forcing.F

C--   Now calculate whether it is time to update the forcing arrays
      CALL GET_PERIODIC_INTERVAL(
     O                  intimeP, intime0, intime1, bWght, aWght,
     I                  externForcingCycle, externForcingPeriod,
     I                  deltaTClock, myTime, myThid )

      bi = myBxLo(myThid)
      bj = myByLo(myThid)
#ifdef ALLOW_DEBUG
      IF ( debugLevel.GE.debLevB ) THEN
        _BEGIN_MASTER(myThid)
        WRITE(standardMessageUnit,'(A,I10,A,4I5,A,2F14.10)')
     &   ' EXTERNAL_FIELDS_LOAD,', myIter,
     &   ' : iP,iLd,i0,i1=', intimeP,loadedRec(bi,bj), intime0,intime1,
     &   ' ; Wght=', bWght, aWght
        _END_MASTER(myThid)
      ENDIF
#endif /* ALLOW_DEBUG */
#ifdef ALLOW_AUTODIFF
C-    assuming that we call S/R EXTERNAL_FIELDS_LOAD at each time-step and
C     with increasing time, this will catch when we need to load new records;
C     But with Adjoint run, this is not always the case => might end-up using
C     the wrong time-records
# ifndef STORE_LOADEDREC_TEST
      IF ( intime0.NE.intimeP .OR. myIter.EQ.nIter0 ) THEN
# else
      IF ( intime1.NE.loadedRec(bi,bj) ) THEN
# endif
#else /* ALLOW_AUTODIFF */
C-    Make no assumption on sequence of calls to EXTERNAL_FIELDS_LOAD ;
C     This is the correct formulation (works in Adjoint run).
C     Unfortunatly, produces many recomputations <== not used until it is fixed
      IF ( intime1.NE.loadedRec(bi,bj) ) THEN
#endif /* ALLOW_AUTODIFF */

C--   If the above condition is met then we need to read in
C     data for the period ahead and the period behind myTime.
        IF ( debugLevel.GE.debLevZero ) THEN
         _BEGIN_MASTER(myThid)
         WRITE(standardMessageUnit,'(A,I10,A,2(2I5,A))')
     &    ' EXTERNAL_FIELDS_LOAD, it=', myIter,
     &    ' : Reading new data, i0,i1=', intime0, intime1,
     &    ' (prev=', intimeP, loadedRec(bi,bj), ' )'
         _END_MASTER(myThid)
        ENDIF

        IF ( zonalWindFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS( zonalWindFile, taux0,
     &                        intime0, myIter, myThid )
         CALL READ_REC_XY_RS( zonalWindFile, taux1,
     &                        intime1, myIter, myThid )
        ENDIF
        IF ( meridWindFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS( meridWindFile, tauy0,
     &                        intime0, myIter, myThid )
         CALL READ_REC_XY_RS( meridWindFile, tauy1,
     &                        intime1, myIter, myThid )
        ENDIF
        IF ( surfQFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS( surfQFile, Qnet0,
     &                        intime0, myIter, myThid )
         CALL READ_REC_XY_RS( surfQFile, Qnet1,
     &                        intime1, myIter, myThid )
        ELSEIF ( surfQnetFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS( surfQnetFile, Qnet0,
     &                        intime0, myIter, myThid )
         CALL READ_REC_XY_RS( surfQnetFile, Qnet1,
     &                        intime1, myIter, myThid )
        ENDIF
        IF ( EmPmRfile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS( EmPmRfile, EmPmR0,
     &                        intime0, myIter, myThid )
         CALL READ_REC_XY_RS( EmPmRfile, EmPmR1,
     &                        intime1, myIter, myThid )
c        IF ( convertEmP2rUnit.EQ.mass2rUnit ) THEN
C-    EmPmR is now (after c59h) expressed in kg/m2/s (fresh water mass flux)
          DO bj = myByLo(myThid), myByHi(myThid)
           DO bi = myBxLo(myThid), myBxHi(myThid)
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
              EmPmR0(i,j,bi,bj) = EmPmR0(i,j,bi,bj)*rhoConstFresh
              EmPmR1(i,j,bi,bj) = EmPmR1(i,j,bi,bj)*rhoConstFresh
             ENDDO
            ENDDO
           ENDDO
          ENDDO
c        ENDIF
        ENDIF
        IF ( saltFluxFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS( saltFluxFile, saltFlux0,
     &                        intime0, myIter, myThid )
         CALL READ_REC_XY_RS( saltFluxFile, saltFlux1,
     &                        intime1, myIter, myThid )
        ENDIF
        IF ( thetaClimFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS( thetaClimFile, SST0,
     &                        intime0, myIter, myThid )
         CALL READ_REC_XY_RS( thetaClimFile, SST1,
     &                        intime1, myIter, myThid )
        ENDIF
        IF ( saltClimFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS( saltClimFile, SSS0,
     &                        intime0, myIter, myThid )
         CALL READ_REC_XY_RS( saltClimFile, SSS1,
     &                        intime1, myIter, myThid )
        ENDIF
#ifdef SHORTWAVE_HEATING
        IF ( surfQswFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS( surfQswFile, Qsw0,
     &                        intime0, myIter, myThid )
         CALL READ_REC_XY_RS( surfQswFile, Qsw1,
     &                        intime1, myIter, myThid )
         IF ( surfQFile .NE. ' '  ) THEN
C-    Qnet is now (after c54) the net Heat Flux (including SW)
          DO bj = myByLo(myThid), myByHi(myThid)
           DO bi = myBxLo(myThid), myBxHi(myThid)
            DO j=1-OLy,sNy+OLy
             DO i=1-OLx,sNx+OLx
              Qnet0(i,j,bi,bj) = Qnet0(i,j,bi,bj) + Qsw0(i,j,bi,bj)
              Qnet1(i,j,bi,bj) = Qnet1(i,j,bi,bj) + Qsw1(i,j,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDIF
#endif
#ifdef ALLOW_GEOTHERMAL_FLUX
        IF ( geothermalFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS( geothermalFile, geothFlux0,
     &                        intime0, myIter, myThid )
         CALL READ_REC_XY_RS( geothermalFile, geothFlux1,
     &                        intime1, myIter, myThid )
         _EXCH_XY_RS( geothFlux0, myThid )
         _EXCH_XY_RS( geothFlux1, myThid )
        ENDIF
#endif
#ifdef ATMOSPHERIC_LOADING
        IF ( pLoadFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS( pLoadFile, pLoad0,
     &                        intime0, myIter, myThid )
         CALL READ_REC_XY_RS( pLoadFile, pLoad1,
     &                        intime1, myIter, myThid )
        ENDIF
#endif

C-    thread synchronisation (barrier) is part of the EXCH S/R calls
        _EXCH_XY_RS(SST0  , myThid )
        _EXCH_XY_RS(SST1  , myThid )
        _EXCH_XY_RS(SSS0  , myThid )
        _EXCH_XY_RS(SSS1  , myThid )
        CALL EXCH_UV_XY_RS(taux0,tauy0,.TRUE.,myThid)
        CALL EXCH_UV_XY_RS(taux1,tauy1,.TRUE.,myThid)
        _EXCH_XY_RS(Qnet0, myThid )
        _EXCH_XY_RS(Qnet1, myThid )
        _EXCH_XY_RS(EmPmR0, myThid )
        _EXCH_XY_RS(EmPmR1, myThid )
        _EXCH_XY_RS(saltFlux0, myThid )
        _EXCH_XY_RS(saltFlux1, myThid )
#ifdef SHORTWAVE_HEATING
        _EXCH_XY_RS(Qsw0, myThid )
        _EXCH_XY_RS(Qsw1, myThid )
#endif
#ifdef ATMOSPHERIC_LOADING
        _EXCH_XY_RS(pLoad0, myThid )
        _EXCH_XY_RS(pLoad1, myThid )
#endif

C-    save newly loaded time-record
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           loadedRec(bi,bj) = intime1
         ENDDO
        ENDDO

C--   end if-block for loading new time-records
      ENDIF

C--   Interpolate fu,fv,Qnet,EmPmR,SST,SSS,Qsw
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        IF ( thetaClimFile .NE. ' '  ) THEN
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            SST(i,j,bi,bj)   = bWght*SST0(i,j,bi,bj)
     &                       + aWght*SST1(i,j,bi,bj)
           ENDDO
          ENDDO
        ENDIF
        IF ( saltClimFile .NE. ' '  ) THEN
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            SSS(i,j,bi,bj)   = bWght*SSS0(i,j,bi,bj)
     &                       + aWght*SSS1(i,j,bi,bj)
           ENDDO
          ENDDO
        ENDIF
        IF ( zonalWindFile .NE. ' '  ) THEN
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            fu(i,j,bi,bj)    = bWght*taux0(i,j,bi,bj)
     &                       + aWght*taux1(i,j,bi,bj)
           ENDDO
          ENDDO
        ENDIF
        IF ( meridWindFile .NE. ' '  ) THEN
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            fv(i,j,bi,bj)    = bWght*tauy0(i,j,bi,bj)
     &                       + aWght*tauy1(i,j,bi,bj)
           ENDDO
          ENDDO
        ENDIF
        IF ( surfQnetFile .NE. ' '
     &     .OR. surfQFile .NE. ' '  ) THEN
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            Qnet(i,j,bi,bj)  = bWght*Qnet0(i,j,bi,bj)
     &                       + aWght*Qnet1(i,j,bi,bj)
           ENDDO
          ENDDO
        ENDIF
        IF ( EmPmRfile .NE. ' '  ) THEN
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            EmPmR(i,j,bi,bj) = bWght*EmPmR0(i,j,bi,bj)
     &                       + aWght*EmPmR1(i,j,bi,bj)
           ENDDO
          ENDDO
        ENDIF
        IF ( saltFluxFile .NE. ' '  ) THEN
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            saltFlux(i,j,bi,bj) = bWght*saltFlux0(i,j,bi,bj)
     &                          + aWght*saltFlux1(i,j,bi,bj)
           ENDDO
          ENDDO
        ENDIF
#ifdef SHORTWAVE_HEATING
        IF ( surfQswFile .NE. ' '  ) THEN
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            Qsw(i,j,bi,bj)   = bWght*Qsw0(i,j,bi,bj)
     &                       + aWght*Qsw1(i,j,bi,bj)
           ENDDO
          ENDDO
        ENDIF
#endif
#ifdef ALLOW_GEOTHERMAL_FLUX
        IF ( geothermalFile .NE. ' '  ) THEN
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            geothermalFlux(i,j,bi,bj) = bWght*geothFlux0(i,j,bi,bj)
     &                                + aWght*geothFlux1(i,j,bi,bj)
           ENDDO
          ENDDO
        ENDIF
#endif
#ifdef ATMOSPHERIC_LOADING
        IF ( pLoadFile .NE. ' '  ) THEN
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            pLoad(i,j,bi,bj) = bWght*pLoad0(i,j,bi,bj)
     &                       + aWght*pLoad1(i,j,bi,bj)
           ENDDO
          ENDDO
        ENDIF
#endif
       ENDDO
      ENDDO

C-- Print for checking:
#ifdef ALLOW_DEBUG
      IF ( debugLevel.GE.debLevC ) THEN
        _BEGIN_MASTER( myThid )
        WRITE(standardMessageUnit,'(A,1P4E12.4)')
     &   ' EXTERNAL_FIELDS_LOAD: (fu0,1),fu,fv=',
     &        taux0(1,sNy,1,1), taux1(1,sNy,1,1),
     &           fu(1,sNy,1,1),    fv(1,sNy,1,1)
        WRITE(standardMessageUnit,'(A,1P4E12.4)')
     &   ' EXTERNAL_FIELDS_LOAD: SST,SSS,Q,E-P=',
     &          SST(1,sNy,1,1),   SSS(1,sNy,1,1),
     &         Qnet(1,sNy,1,1), EmPmR(1,sNy,1,1)
        _END_MASTER( myThid )
      ENDIF
#endif /* ALLOW_DEBUG */

C endif for periodicForcing
      ENDIF

CHL   New block for SGS forcing
      IF ( sgsExternalForcing ) THEN

C--   Now calculate whether it is time to update the forcing arrays
      CALL GET_PERIODIC_INTERVAL(
     O                  sgsIntimeP, sgsIntime0, sgsIntime1,
     O                  sgsBWght, sgsAWght,
     I                  externForcingCycle, externForcingPeriod,
     I                  deltaTClock,
     I                  myTime + 0.5 _d 0 * externForcingPeriod,
     I                  myThid )

      sgsBWght = 1. _d 0 - sgsBWght
      sgsAWght = 1. _d 0 - sgsAWght

#ifdef ALLOW_DEBUG
      WRITE(standardMessageUnit,'(A,I10,A,3I5,A,2F14.10)')
     &   ' EXTERNAL_FIELDS_LOAD,', myIter,
     &   ' : iP,i0,i1=', sgsIntimeP,sgsIntime0,sgsIntime1,
     &   ' ; Wght=', sgsBWght,sgsAWght
#endif /* ALLOW_DEBUG */

      IF ( sgsIntime0 .NE. sgsIntimeP ) THEN
        PRINT *, "reading the forcing files!"
        IF ( sgsThetaFile .NE. ' ' ) THEN
         CALL READ_REC_XYZ_RS( sgsThetaFile, sgsT0,
     &                         sgsIntime0, myIter, myThid )
         CALL READ_REC_XYZ_RS( sgsThetaFile, sgsT1,
     &                         sgsIntime1, myIter, myThid )
        ENDIF

        IF ( sgsSaltFile .NE. ' ' ) THEN
          CALL READ_REC_XYZ_RS(
     &         sgsSaltFile, sgsS0,
     &         sgsIntime0, myIter, myThid )
          CALL READ_REC_XYZ_RS(
     &         sgsSaltFile, sgsS1,
     &         sgsIntime1, myIter, myThid )
        ENDIF

C-    thread synchronisation (barrier) is part of the EXCH S/R calls
        _EXCH_XYZ_RS( sgsT0, myThid )
        _EXCH_XYZ_RS( sgsT1, myThid )
        _EXCH_XYZ_RS( sgsS0, myThid )
        _EXCH_XYZ_RS( sgsS1, myThid )

C--   end if-block for loading new time-records
      ENDIF

C--   Interpolate fu,fv,Qnet,EmPmR,SST,SSS,Qsw
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        IF ( sgsThetaFile .NE. ' ' ) THEN
          DO k=1,Nr
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             sgsT(i,j,k,bi,bj) =
     &          sgsBWght*sgsT0(i,j,k,bi,bj)
     &        + sgsAWght*sgsT1(i,j,k,bi,bj)
            ENDDO
           ENDDO
          ENDDO
        ENDIF
        IF ( sgsSaltFile .NE. ' ' ) THEN
          DO k=1,Nr
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             sgsS(i,j,k,bi,bj) =
     &          sgsBWght*sgsS0(i,j,k,bi,bj)
     &        + sgsAWght*sgsS1(i,j,k,bi,bj)
            ENDDO
           ENDDO
          ENDDO
        ENDIF
       ENDDO
      ENDDO

#ifdef ALLOW_DEBUG
      IF ( debugLevel.GE.debLevC ) THEN
        _BEGIN_MASTER( myThid )

        IF ( sgsThetaFile .NE. ' ' ) THEN
          WRITE(standardMessageUnit,'(A,1P6E12.4)')
     &      ' EXTERNAL_FIELDS_LOAD: SGS T=',
     &      sgsT0(10,10,10,1,1),
     &      sgsT1(10,10,10,1,1),
     &      sgsT (10,10,10,1,1),
     &      sgsBWght, sgsAWght, 0. _d 0
        ENDIF
        IF ( sgsSaltFile .NE. ' ' ) THEN
          WRITE(standardMessageUnit,'(A,1P6E12.4)')
     &      ' EXTERNAL_FIELDS_LOAD: SGS S=',
     &      sgsS0(10,10,10,1,1),
     &      sgsS1(10,10,10,1,1),
     &      sgsS (10,10,10,1,1),
     &      sgsBWght, sgsAWght, 0. _d 0
        ENDIF

        _END_MASTER( myThid )
      ENDIF
#endif /* ALLOW_DEBUG */

C endif for sgsExternalForcing
      ENDIF

#endif /* EXCLUDE_FFIELDS_LOAD */

      RETURN
      END
