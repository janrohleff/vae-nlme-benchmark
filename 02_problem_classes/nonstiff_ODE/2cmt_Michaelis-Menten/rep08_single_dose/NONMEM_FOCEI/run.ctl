$PROBLEM 2-cmt Michaelis-Menten | pk2mm rep08 | FOCEI
$ABBR DERIV2=NO
$INPUT ID TIME AMT DV MDV EVID CMT
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=tolswitch.f90
$MODEL COMP=(CENTRAL,DEFDOSE,DEFOBS) COMP=(PERIPH)

$PK
; individual parameters: theta_i = exp(MU + ETA), MU = log(typical value) [+ covariate effects];
; parameters without IIV have an ETA with OMEGA 0 FIX (MU referencing for FOCE FAST)
MU_1 = LOG(THETA(1))
V1 = EXP(MU_1 + ETA(1))
MU_2 = LOG(THETA(2))
Q = EXP(MU_2 + ETA(2))
MU_3 = LOG(THETA(3))
V2 = EXP(MU_3 + ETA(3))
MU_4 = LOG(THETA(4))
VM = EXP(MU_4 + ETA(4))
MU_5 = LOG(THETA(5))
KM = EXP(MU_5 + ETA(5))

$DES
CC = A(1)/V1
DADT(1) = -VM*CC/(KM + CC) - Q/V1*A(1) + Q/V2*A(2)
DADT(2) = Q/V1*A(1) - Q/V2*A(2)

$ERROR
IPRED = A(1)/V1
Y = IPRED + IPRED*EPS(1)

$THETA
(0, 14.0)   ; 1 V1_pop
(0, 3.0)   ; 2 Q_pop
(0, 21.0)   ; 3 V2_pop
(0, 14.0)   ; 4 Vm_pop
(0, 3.0)   ; 5 Km_pop

$OMEGA
0.25   ; ETA(1) V1  (omega = 0.5)
0 FIX   ; ETA(2) Q  (no IIV)
0.25   ; ETA(3) V2  (omega = 0.5)
0.25   ; ETA(4) Vm  (omega = 0.5)
0.25   ; ETA(5) Km  (omega = 0.5)

$SIGMA
0.09   ; EPS(1) Conc proportional  (b = 0.3)

; -2LL: importance sampling (EONLY=1), proposal t with 5 df; first IMP adapts the proposal,
; the second evaluates with 10000 samples per individual (as Monolix). Reported value:
; 'OBJECTIVE FUNCTION VALUE WITH CONSTANT' of the last step.
$EST METHOD=COND INTERACTION FAST MAXEVAL=9999 NSIG=2 SIGL=6 PRINT=5 NOABORT
     FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=5 ISAMPLE=1000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=23456 NOCOV=1 FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=1 ISAMPLE=10000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=34567 NOCOV=1 FNLETA=0
