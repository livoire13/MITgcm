CBOP
C    !ROUTINE: ECCO.h
C    !INTERFACE:
C #include ECCO.h

C    !DESCRIPTION: \bv
C     ==================================================================
C     ECCO.h
C     o Main header file for the ECCO pkg.
C       started: Christian Eckert eckert@mit.edu  23-Feb-1999
C     ==================================================================
C     \ev
CEOP

C     Version of the release and versions of used packages:
C     =====================================================
C
C     eccoVersion                - ecco release version.
C     usesCalendarVersion        - version of the calendar that has to
C                                  be used.
C     usesExternalForcingVersion - version of the external forcing that
C                                  has to be used.
C     usesAdjointSupportVersion  - version of the adjoint support routines
C                                  that have to be used.
C     usesOptimizationVersion    - version of the oof-line optimization
C                                  that has to be used.

      CHARACTER*(5) eccoVersion
      CHARACTER*(5) eccoUsesCalVersion
      CHARACTER*(5) eccoUsesExfVersion
      CHARACTER*(5) eccoUsesAdsuppVersion
      CHARACTER*(5) eccoUsesOptimVersion

      PARAMETER(    eccoVersion           = '0.1.0' )
      PARAMETER(    eccoUsesCalVersion    = '0.1.4' )
      PARAMETER(    eccoUsesExfVersion    = '0.1.1' )
      PARAMETER(    eccoUsesAdsuppVersion = '0.1.0' )
      PARAMETER(    eccoUsesOptimVersion  = '2.1.0' )

C     Experiment name:
C     ================

      COMMON /ECCO_C/
     &                expId
      CHARACTER*(10)  expId

C     Integration information:
C     ========================
C
C     nyears - number of calendar years that are affected by the
C              current integration.

      COMMON /ECCO_I/
     &                nyears, nmonths, ndays, numsteps,
     &                eccoiter, ecco_prevcall
      INTEGER nyears
      INTEGER nmonths
      INTEGER ndays
      INTEGER numsteps
      INTEGER eccoiter
      INTEGER ecco_prevcall

C     Averaging counters:
C     ===================
C
C     sum1day :: counter for the daily averaging
C     sum1mon :: counter for the monthly averaging
C     dayrec  :: number of averaged surface pressure records.
C     monrec  :: number of averaged theta and salinity records.

      COMMON /AVERAGE_I/
     &                   sum1day,sum1mon,sum1year,
     &                   dayrec,monrec,yearrec
      INTEGER sum1day
      INTEGER sum1mon
      INTEGER sum1year
      INTEGER dayrec
      INTEGER monrec
      INTEGER yearrec

C     Flags used in the model-data comparison:
C     ========================================
C
C     using_ers :: flag that indicates the use of ERS data

      COMMON /ECCO_COST_DATA_FLAGS/
     &                         using_mdt,
     &                         using_tpj,
     &                         using_topex,
     &                         using_ers,
     &                         using_gfo,
     &                         using_cost_altim,
     &                         using_cost_sst,
     &                         using_cost_seaice,
     &                         using_cost_transp
      LOGICAL using_mdt
      LOGICAL using_tpj
      LOGICAL using_topex
      LOGICAL using_ers
      LOGICAL using_gfo
      LOGICAL using_cost_altim
      LOGICAL using_cost_sst
      LOGICAL using_cost_seaice
      LOGICAL using_cost_transp

C     Record counters relevant for the cost function evaluation.
C     ==========================================================
C
C     nyearsrec:: number of yearly records that will be generated by
C                 the current model integration.
C     nmonsrec :: number of monthly records that will be generated by
C                 the current model integration.
C     ndaysrec :: number of  daily  records that will be generated by
C                 the current model integration.

      COMMON /ECCO_I/
     &                nyearsrec,
     &                nmonsrec,
     &                ndaysrec
      INTEGER nyearsrec
      INTEGER nmonsrec
      INTEGER ndaysrec

      COMMON /ECCO_R/
     &                    m_eta, m_UE, m_VN, m_bp,
     &                    trVol, trHeat, trSalt,
     &                    frame
      _RL m_eta    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,   nSx,nSy)
      _RL m_UE     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL m_VN     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL m_bp     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,   nSx,nSy)
      _RL trVol    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL trHeat   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL trSalt   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL frame    (1-OLx:sNx+OLx,1-OLy:sNy+OLy           )

#ifdef ALLOW_PSBAR_STERIC
      COMMON /ECCO_R2/
     &                    sterGloH,
     &                    VOLsumGlob_0, VOLsumGlob,
     &                    RHOsumGlob_0, RHOsumGlob
      _RL sterGloH
      _RL VOLsumGlob_0, VOLsumGlob, RHOsumGlob_0, RHOsumGlob

#endif
#ifdef ATMOSPHERIC_LOADING
#ifdef ALLOW_IB_CORR
      COMMON /ECCO_R3/
     &                    m_bp_nopabar,
     &                    m_eta_dyn
      _RL m_bp_nopabar(1-OLx:sNx+OLx,1-OLy:sNy+OLy,   nSx,nSy)
      _RL m_eta_dyn   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,   nSx,nSy)

#endif
#endif
#ifndef ECCO_VARIABLE_AREAVOLGLOB
      COMMON /ECCO_R4/ eccoVol_0
      _RL eccoVol_0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
#endif

C     Two runtime parameters related to outputting sterGloH
C     ecco_output_sterGloH :: output sterGloH at each time step if true
C     ecco_keepTSeriesOutp_open :: keep the sterGloH file open if true
C     eccoWriteCostFunction :: internal flag to suppress writing
C                       costfunction multiple times (e.g. with grdchk-pkg)
      COMMON /ECCO_L/
     &     ecco_output_sterGloH, ecco_keepTSeriesOutp_open,
     &     eccoWriteCostFunction
      LOGICAL ecco_output_sterGloH, ecco_keepTSeriesOutp_open
      LOGICAL eccoWriteCostFunction

C     file precision
      COMMON /PREC_TYPE_COST/
     &                        cost_iprec
      INTEGER cost_iprec

C     empty pre/post-processing :
C     =========================
      COMMON /ECCO_NOGENCOST_C/
     &       no_preproc,
     &       no_preproc_c,
     &       no_posproc,
     &       no_posproc_c,
     &       clim_preproc,
     &       anom_preproc
      CHARACTER*(MAX_LEN_FNAM) no_preproc(NGENPPROC)
      CHARACTER*(MAX_LEN_FNAM) no_preproc_c(NGENPPROC)
      CHARACTER*(MAX_LEN_FNAM) no_posproc(NGENPPROC)
      CHARACTER*(MAX_LEN_FNAM) no_posproc_c(NGENPPROC)
      CHARACTER*(MAX_LEN_FNAM) clim_preproc(NGENPPROC)
      CHARACTER*(MAX_LEN_FNAM) anom_preproc(NGENPPROC)

      COMMON /ECCO_NOGENCOST_R/
     &       no_preproc_r, no_posproc_r
      _RL no_preproc_r(NGENPPROC)
      _RL no_posproc_r(NGENPPROC)

      COMMON /ECCO_NOGENCOST_I/
     &       no_preproc_i, no_posproc_i, clim_preproc_i
      INTEGER no_preproc_i(NGENPPROC)
      INTEGER no_posproc_i(NGENPPROC)
      INTEGER clim_preproc_i(NGENPPROC)

C     gencost common blocs:
C     =====================

#ifdef ALLOW_GENCOST_CONTRIBUTION
C     objf_gencost :: gencost user defined contribution
      COMMON /ECCO_GENCOST_CTRL/
     &       gencost_dummy
      _RL  gencost_dummy(NGENCOST)

      COMMON /ECCO_GENCOST_R_1/
     &       objf_gencost, num_gencost, mult_gencost, gencost_storefld,
     &       gencost_barfld, gencost_modfld, gencost_weight,
     &       gencost_mskCsurf, gencost_mskWsurf, gencost_mskSsurf,
     &       gencost_mskVertical,
#ifdef ALLOW_GENCOST3D
     &       gencost_bar3d, gencost_mod3d, gencost_wei3d,
     &       gencost_mskC, gencost_mskW, gencost_mskS,
#endif
     &       gencost_refPressure, gencost_sigmaLow, gencost_sigmaHigh,
     &       gencost_tanhScale,
     &       gencost_spmin, gencost_spmax, gencost_spzero,
     &       gencost_period, gencost_preproc_r, gencost_posproc_r,
     &       gencost_wei1d, gencost_1ddata

      _RL  objf_gencost(nSx,nSy,NGENCOST)
      _RL  num_gencost(nSx,nSy,NGENCOST)
      _RL  mult_gencost(NGENCOST)
      _RL  gencost_spmin(NGENCOST)
      _RL  gencost_spmax(NGENCOST)
      _RL  gencost_spzero(NGENCOST)
      _RL  gencost_period(NGENCOST)
      _RL  gencost_storefld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       nSx,nSy,NGENCOST)
      _RL  gencost_barfld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       nSx,nSy,NGENCOST)
      _RL  gencost_modfld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       nSx,nSy,NGENCOST)
      _RL  gencost_weight(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       nSx,nSy,NGENCOST)
      _RL  gencost_mskCsurf(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       nSx,nSy,NGENCOST)
      _RL  gencost_mskWsurf(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       nSx,nSy,NGENCOST)
      _RL  gencost_mskSsurf(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       nSx,nSy,NGENCOST)
      _RL  gencost_mskVertical(Nr,NGENCOST)
      _RL  gencost_sigmaLow(NGENCOST)
      _RL  gencost_sigmaHigh(NGENCOST)
      _RL  gencost_refPressure(NGENCOST)
      _RL  gencost_tanhScale(NGENCOST)
#ifdef ALLOW_GENCOST3D
      _RL  gencost_bar3d(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       Nr,nSx,nSy,NGENCOST3D)
      _RL  gencost_mod3d(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       Nr,nSx,nSy,NGENCOST3D)
      _RL  gencost_wei3d(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       Nr,nSx,nSy,NGENCOST3D)
      _RL  gencost_mskC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       Nr,nSx,nSy,NGENCOST3D)
      _RL  gencost_mskW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       Nr,nSx,nSy,NGENCOST3D)
      _RL  gencost_mskS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &       Nr,nSx,nSy,NGENCOST3D)
#endif
      _RL gencost_preproc_r(NGENPPROC,NGENCOST)
      _RL gencost_posproc_r(NGENPPROC,NGENCOST)
      _RL gencost_wei1d(NGENCOST)
      _RL gencost_1ddata(N1DDATA, NGENCOST)

C     gencost_kLev_select :: selected vertical level of a 3D model field to use
C                            to compute a cost term contribution (default = 1).
C                            Note: currently only used in drifter velocity cost

      COMMON /ECCO_GENCOST_I_1/
     &       gencost_nrec, gencost_nrecperiod,
     &       gencost_flag, gencost_outputlevel,
     &       gencost_startdate1, gencost_startdate2,
     &       gencost_enddate1, gencost_enddate2,
     &       gencost_startdate, gencost_enddate,
     &       gencost_pointer3d, gencost_smooth2Ddiffnbt,
     &       gencost_preproc_i, gencost_posproc_i,
     &       gencost_msk_pointer3d, gencost_itracer,
     &       gencost_kLev_select

      INTEGER gencost_nrec(NGENCOST)
      INTEGER gencost_nrecperiod(NGENCOST)
      INTEGER gencost_flag(NGENCOST)
      INTEGER gencost_outputlevel(NGENCOST)
      INTEGER gencost_startdate1(NGENCOST)
      INTEGER gencost_startdate2(NGENCOST)
      INTEGER gencost_startdate(4,NGENCOST)
      INTEGER gencost_enddate1(NGENCOST)
      INTEGER gencost_enddate2(NGENCOST)
      INTEGER gencost_enddate(4,NGENCOST)
      INTEGER gencost_pointer3d(NGENCOST)
      INTEGER  gencost_smooth2Ddiffnbt(NGENCOST)
      INTEGER gencost_preproc_i(NGENPPROC,NGENCOST)
      INTEGER gencost_posproc_i(NGENPPROC,NGENCOST)
      INTEGER gencost_msk_pointer3d(NGENCOST)
      INTEGER gencost_itracer(NGENCOST)
      INTEGER gencost_kLev_select(NGENCOST)

      COMMON /ECCO_GENCOST_L_1/
     &       gencost_timevaryweight, gencost_barskip,
     &       using_gencost, gencost_is3d, gencost_msk_is3d,
     &       gencost_is1d, gencost_useDensityMask
      LOGICAL using_gencost(NGENCOST)
      LOGICAL gencost_is3d(NGENCOST)
      LOGICAL gencost_is1d(NGENCOST)
      LOGICAL gencost_msk_is3d(NGENCOST)
      LOGICAL gencost_useDensityMask(NGENCOST)
      LOGICAL gencost_timevaryweight(NGENCOST)
      LOGICAL gencost_barskip(NGENCOST)

      COMMON /ECCO_GENCOST_C/
     &       gencost_name,
     &       gencost_scalefile,
     &       gencost_errfile,
     &       gencost_datafile,
     &       gencost_barfile,
     &       gencost_avgperiod,
     &       gencost_preproc,
     &       gencost_preproc_c,
     &       gencost_posproc,
     &       gencost_posproc_c,
     &       gencost_mask
      CHARACTER*(MAX_LEN_FNAM) gencost_name(NGENCOST)
      CHARACTER*(MAX_LEN_FNAM) gencost_scalefile(NGENCOST)
      CHARACTER*(MAX_LEN_FNAM) gencost_errfile(NGENCOST)
      CHARACTER*(MAX_LEN_FNAM) gencost_datafile(NGENCOST)
      CHARACTER*(MAX_LEN_FNAM) gencost_barfile(NGENCOST)
      CHARACTER*(5)            gencost_avgperiod(NGENCOST)
      CHARACTER*(MAX_LEN_FNAM) gencost_preproc(NGENPPROC,NGENCOST)
      CHARACTER*(MAX_LEN_FNAM) gencost_posproc_c(NGENPPROC,NGENCOST)
      CHARACTER*(MAX_LEN_FNAM) gencost_posproc(NGENPPROC,NGENCOST)
      CHARACTER*(MAX_LEN_FNAM) gencost_preproc_c(NGENPPROC,NGENCOST)
      CHARACTER*(MAX_LEN_FNAM) gencost_mask(NGENCOST)

#endif /* ALLOW_GENCOST_CONTRIBUTION */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
