$PROBLEM Effect compartment | effcmt N020_rep04 | FOCEI
$ABBR DERIV2=NO
$INPUT ID TIME AMT DV MDV EVID CMT DVID
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=tolswitch.f90
$MODEL COMP=(DEPOT,DEFDOSE) COMP=(CENTRAL,DEFOBS) COMP=(EFFECT)

$PK
; individual parameters: theta_i = exp(MU + ETA), MU = log(typical value) [+ covariate effects];
; parameters without IIV have an ETA with OMEGA 0 FIX (MU referencing for FOCE FAST)
MU_1 = LOG(THETA(1))
KA = EXP(MU_1 + ETA(1))
MU_2 = LOG(THETA(2))
V = EXP(MU_2 + ETA(2))
MU_3 = LOG(THETA(3))
CL = EXP(MU_3 + ETA(3))
MU_4 = LOG(THETA(4))
KE0 = EXP(MU_4 + ETA(4))
MU_5 = LOG(THETA(5))
EBASE = EXP(MU_5 + ETA(5))
MU_6 = LOG(THETA(6))
EMAX = EXP(MU_6 + ETA(6))
MU_7 = LOG(THETA(7))
EC50 = EXP(MU_7 + ETA(7))

$DES
CC = A(2)/V
DADT(1) = -KA*A(1)
DADT(2) = KA*A(1) - CL/V*A(2)
DADT(3) = KE0*(CC - A(3))

$ERROR
IPRED = A(2)/V
IF (DVID.EQ.2) IPRED = EBASE + EMAX*A(3)/(EC50 + A(3))
Y = IPRED + IPRED*EPS(1)
IF (DVID.EQ.2) Y = IPRED + EPS(2)

$THETA
(0, 0.6)   ; 1 ka_pop
(0, 70.0)   ; 2 V_pop
(0, 2.8)   ; 3 Cl_pop
(0, 0.45)   ; 4 ke0_pop
(0, 12.0)   ; 5 E0_pop
(0, 35.0)   ; 6 Emax_pop
(0, 3.0)   ; 7 EC50_pop

$OMEGA
0.25   ; ETA(1) ka  (omega = 0.5)
0.25   ; ETA(2) V  (omega = 0.5)
0.25   ; ETA(3) Cl  (omega = 0.5)
0.25   ; ETA(4) ke0  (omega = 0.5)
0.25   ; ETA(5) E0  (omega = 0.5)
0.25   ; ETA(6) Emax  (omega = 0.5)
0.25   ; ETA(7) EC50  (omega = 0.5)

$SIGMA
0.09   ; EPS(1) Conc proportional  (b1 = 0.3)
25.0   ; EPS(2) Effect additive  (a2 = 5)

; -2LL: importance sampling (EONLY=1), proposal t with 5 df; first IMP adapts the proposal,
; the second evaluates with 10000 samples per individual (as Monolix). Reported value:
; 'OBJECTIVE FUNCTION VALUE WITH CONSTANT' of the last step.
$EST METHOD=COND INTERACTION FAST MAXEVAL=9999 NSIG=2 SIGL=6 PRINT=5 NOABORT
     FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=5 ISAMPLE=1000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=23456 NOCOV=1 FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=1 ISAMPLE=10000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=34567 NOCOV=1 FNLETA=0
