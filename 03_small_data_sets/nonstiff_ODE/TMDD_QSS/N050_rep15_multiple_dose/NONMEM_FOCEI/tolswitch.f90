      SUBROUTINE TOL(NRD,ANRD,NRDC,ANRDC)
! ODE-Toleranz je Schätzmethode (Hilfe $TOL, TOL-Routine; NM_STEP/IBMETHOD aus NONMEMs Modulen):
! Importance Sampling (IBMETHOD = EST_IMP, −2LL) rtol 1e-3 / atol 1e-6 wie Monolix und VAE Studio,
! alles andere (FOCE-Schätzung) rtol 1e-6 / atol 1e-6 nach NONMEM-Handbuch.
      USE EST_DEFS, ONLY: EST_IMP
      USE NM_BAYES_INT, ONLY: NM_STEP, EST_STEP
      USE NMBAYES_INT, ONLY: IBMETHOD
      IMPLICIT NONE
      INTEGER :: NRD(0:*), ANRD(0:*), NRDC(0:*), ANRDC(0:*)
      IF (NM_STEP == EST_STEP .AND. IBMETHOD == EST_IMP) THEN
        NRD(1) = 3
        ANRD(1) = 6
      ELSE
        NRD(1) = 6
        ANRD(1) = 6
      ENDIF
      RETURN
      END
