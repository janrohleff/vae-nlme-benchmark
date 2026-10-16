$PROBLEM Daridorexant | daridorexant rep01 | FOCEI
$ABBR DERIV2=NO
$INPUT ID TIME AMT DV MDV EVID CMT LDOSE FASTED FED MORNING LLBW LFM LAGE LALP
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=9 ATOL=12
$MODEL COMP=(DEPOT,DEFDOSE) COMP=(CENTRAL,DEFOBS) COMP=(PERIPH)

$PK
; individual parameters: theta_i = exp(MU + ETA), MU = log(typical value) [+ covariate effects];
; parameters without IIV have an ETA with OMEGA 0 FIX (MU referencing for FOCE FAST)
MU_1 = LOG(THETA(1)/(1 - THETA(1))) + THETA(9)*LDOSE
FBIO = 1/(1 + EXP(-(MU_1 + ETA(1))))
MU_2 = LOG(THETA(2)) + THETA(10)*FASTED + THETA(11)*FED
TLAG = EXP(MU_2 + ETA(2))
MU_3 = LOG(THETA(3)) + THETA(12)*FASTED + THETA(13)*FED + THETA(14)*MORNING
KA = EXP(MU_3 + ETA(3))
MU_4 = LOG(THETA(4)) + THETA(15)*LLBW
VC = EXP(MU_4 + ETA(4))
MU_5 = LOG(THETA(5)) + THETA(16)*LFM
VP = EXP(MU_5 + ETA(5))
MU_6 = LOG(THETA(6)) + THETA(17)*LFM
Q = EXP(MU_6 + ETA(6))
MU_7 = LOG(THETA(7))
VM = EXP(MU_7 + ETA(7))
MU_8 = LOG(THETA(8)) + THETA(18)*LLBW + THETA(19)*LAGE + THETA(20)*LALP
KM = EXP(MU_8 + ETA(8))
F1 = FBIO
ALAG1 = TLAG

$DES
CC = A(2)/VC
DADT(1) = -KA*A(1)
DADT(2) = KA*A(1) - VM*CC/(KM + CC) - Q/VC*A(2) + Q/VP*A(3)
DADT(3) = Q/VC*A(2) - Q/VP*A(3)

$ERROR
IPRED = 1000*A(2)/VC
W = THETA(21) + THETA(22)*IPRED
Y = IPRED + W*EPS(1)

$THETA
(0, 0.49029543, 1)   ; 1 F_pop
(0, 0.28444636)   ; 2 tlag_pop
(0, 0.63113737)   ; 3 ka_pop
(0, 20.464648)   ; 4 Vc_pop
(0, 9.5717287)   ; 5 Vp_pop
(0, 5.0075225)   ; 6 Q_pop
(0, 4.8546334)   ; 7 Vm_pop
(0, 3.5382748)   ; 8 Km_pop
0.1   ; 9 beta_F_COV_LDOSE
0.1   ; 10 beta_tlag_COV_FASTED
0.1   ; 11 beta_tlag_COV_FED
0.1   ; 12 beta_ka_COV_FASTED
0.1   ; 13 beta_ka_COV_FED
0.1   ; 14 beta_ka_COV_MORNING
0.1   ; 15 beta_Vc_COV_LLBW
0.1   ; 16 beta_Vp_COV_LFM
0.1   ; 17 beta_Q_COV_LFM
0.1   ; 18 beta_Km_COV_LLBW
0.1   ; 19 beta_Km_COV_LAGE
0.1   ; 20 beta_Km_COV_LALP
(0, 10.0)   ; 21 a
(0, 0.3)   ; 22 b

$OMEGA
0.25   ; ETA(1) F  (omega = 0.5)
0.25   ; ETA(2) tlag  (omega = 0.5)
0.25   ; ETA(3) ka  (omega = 0.5)
0.25   ; ETA(4) Vc  (omega = 0.5)
0.25   ; ETA(5) Vp  (omega = 0.5)
0.25   ; ETA(6) Q  (omega = 0.5)
0.25   ; ETA(7) Vm  (omega = 0.5)
0.25   ; ETA(8) Km  (omega = 0.5)

$SIGMA
1 FIX   ; EPS(1) Conc, sd = a + b*IPRED via THETA

; -2LL: importance sampling (EONLY=1), proposal t with 5 df; first IMP adapts the proposal,
; the second evaluates with 10000 samples per individual (as Monolix). Reported value:
; 'OBJECTIVE FUNCTION VALUE WITH CONSTANT' of the last step.
$EST METHOD=COND INTERACTION FAST MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT
     FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=5 ISAMPLE=1000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=9 PRINT=1 SEED=23456 NOCOV=1 FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=1 ISAMPLE=10000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=9 PRINT=1 SEED=34567 NOCOV=1 FNLETA=0
