$PROBLEM 1-cmt Michaelis-Menten | pkmm rep17 | FOCEI
$ABBR DERIV2=NO
$INPUT ID TIME AMT DV MDV EVID CMT
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=3 ATOL=6
$MODEL COMP=(DEPOT,DEFDOSE) COMP=(CENTRAL,DEFOBS)

$PK
; individual parameters: theta_i = exp(MU + ETA), MU = log(typical value) [+ covariate effects];
; parameters without IIV have an ETA with OMEGA 0 FIX (MU referencing for FOCE FAST)
MU_1 = LOG(THETA(1))
KA = EXP(MU_1 + ETA(1))
MU_2 = LOG(THETA(2))
V = EXP(MU_2 + ETA(2))
MU_3 = LOG(THETA(3))
VM = EXP(MU_3 + ETA(3))
MU_4 = LOG(THETA(4))
KM = EXP(MU_4 + ETA(4))

$DES
CC = A(2)/V
DADT(1) = -KA*A(1)
DADT(2) = KA*A(1) - VM*CC/(KM + CC)

$ERROR
IPRED = A(2)/V
Y = IPRED + IPRED*EPS(1)

$THETA
(0, 0.6)   ; 1 ka_pop
(0, 70.0)   ; 2 V_pop
(0, 7.0)   ; 3 Vm_pop
(0, 3.0)   ; 4 Km_pop

$OMEGA
0.25   ; ETA(1) ka  (omega = 0.5)
0.25   ; ETA(2) V  (omega = 0.5)
0.25   ; ETA(3) Vm  (omega = 0.5)
0.25   ; ETA(4) Km  (omega = 0.5)

$SIGMA
0.09   ; EPS(1) Conc proportional  (b = 0.3)

; -2LL: importance sampling (EONLY=1), proposal t with 5 df; first IMP adapts the proposal,
; the second evaluates with 10000 samples per individual (as Monolix). Reported value:
; 'OBJECTIVE FUNCTION VALUE WITH CONSTANT' of the last step.
$EST METHOD=COND INTERACTION FAST MAXEVAL=9999 NSIG=1 SIGL=3 PRINT=5 NOABORT
     FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=5 ISAMPLE=1000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=23456 NOCOV=1 FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=1 ISAMPLE=10000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=34567 NOCOV=1 FNLETA=0
