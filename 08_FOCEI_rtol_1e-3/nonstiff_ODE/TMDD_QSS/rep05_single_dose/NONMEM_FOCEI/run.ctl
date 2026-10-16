$PROBLEM TMDD (QSS) | tmdd rep05 | FOCEI
$ABBR DERIV2=NO
$INPUT ID TIME AMT DV MDV EVID CMT DVID
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=3 ATOL=6
$MODEL COMP=(ATOT,DEFDOSE,DEFOBS) COMP=(PERIPH) COMP=(RTOT)

$PK
; individual parameters: theta_i = exp(MU + ETA), MU = log(typical value) [+ covariate effects];
; parameters without IIV have an ETA with OMEGA 0 FIX (MU referencing for FOCE FAST)
MU_1 = LOG(THETA(1))
CL = EXP(MU_1 + ETA(1))
MU_2 = LOG(THETA(2))
V1 = EXP(MU_2 + ETA(2))
MU_3 = LOG(THETA(3))
Q = EXP(MU_3 + ETA(3))
MU_4 = LOG(THETA(4))
V2 = EXP(MU_4 + ETA(4))
MU_5 = LOG(THETA(5))
RBASE = EXP(MU_5 + ETA(5))
MU_6 = LOG(THETA(6))
KDEG = EXP(MU_6 + ETA(6))
MU_7 = LOG(THETA(7))
KINT = EXP(MU_7 + ETA(7))
MU_8 = LOG(THETA(8))
KSS = EXP(MU_8 + ETA(8))
A_0(3) = RBASE

$DES
CTOT = A(1)/V1
DQ = CTOT - A(3) - KSS
CF = 0.5*(DQ + SQRT(DQ*DQ + 4*KSS*CTOT))
OCC = A(3)*CF/(KSS + CF)
DADT(1) = -CL*CF - Q*CF + Q*A(2)/V2 - KINT*V1*OCC
DADT(2) = Q*CF - Q*A(2)/V2
DADT(3) = KDEG*RBASE - KDEG*A(3) - (KINT - KDEG)*OCC

$ERROR
CTOTE = A(1)/V1
DQE = CTOTE - A(3) - KSS
CFE = 0.5*(DQE + SQRT(DQE*DQE + 4*KSS*CTOTE))
IPRED = CFE
IF (DVID.EQ.2) IPRED = A(3)
Y = IPRED + IPRED*EPS(1)
IF (DVID.EQ.2) Y = IPRED + IPRED*EPS(2)

$THETA
(0, 0.28)   ; 1 Cl_pop
(0, 2.1)   ; 2 V1_pop
(0, 0.75)   ; 3 Q_pop
(0, 1.75)   ; 4 V2_pop
(0, 14.0)   ; 5 R0_pop
(0, 0.35)   ; 6 kdeg_pop
(0, 0.075)   ; 7 kint_pop
(0, 0.6)   ; 8 KSS_pop

$OMEGA
0.25   ; ETA(1) Cl  (omega = 0.5)
0.25   ; ETA(2) V1  (omega = 0.5)
0 FIX   ; ETA(3) Q  (no IIV)
0.25   ; ETA(4) V2  (omega = 0.5)
0.25   ; ETA(5) R0  (omega = 0.5)
0 FIX   ; ETA(6) kdeg  (no IIV)
0 FIX   ; ETA(7) kint  (no IIV)
0 FIX   ; ETA(8) KSS  (no IIV)

$SIGMA
0.09   ; EPS(1) Free proportional  (b1 = 0.3)
0.09   ; EPS(2) Target proportional  (b2 = 0.3)

; -2LL: importance sampling (EONLY=1), proposal t with 5 df; first IMP adapts the proposal,
; the second evaluates with 10000 samples per individual (as Monolix). Reported value:
; 'OBJECTIVE FUNCTION VALUE WITH CONSTANT' of the last step.
$EST METHOD=COND INTERACTION FAST MAXEVAL=9999 NSIG=1 SIGL=3 PRINT=5 NOABORT
     FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=5 ISAMPLE=1000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=23456 NOCOV=1 FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=1 ISAMPLE=10000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=34567 NOCOV=1 FNLETA=0
