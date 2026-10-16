$PROBLEM Schoemaker bolus_1cpt_VVmKm | sch:bolus_1cpt_VVmKm runN010 | FOCEI
$ABBR DERIV2=NO
$INPUT ID TIME AMT DV MDV EVID CMT
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=tolswitch.f90
$MODEL COMP=(CENTRAL,DEFDOSE,DEFOBS)

$PK
; individual parameters: theta_i = exp(MU + ETA), MU = log(typical value) [+ covariate effects];
; parameters without IIV have an ETA with OMEGA 0 FIX (MU referencing for FOCE FAST)
MU_1 = LOG(THETA(1))
VM = EXP(MU_1 + ETA(1))
MU_2 = LOG(THETA(2))
KM = EXP(MU_2 + ETA(2))
MU_3 = LOG(THETA(3))
V = EXP(MU_3 + ETA(3))

$DES
CC = A(1)/V
DADT(1) = -VM*CC/(KM + CC)

$ERROR
IPRED = A(1)/V
Y = IPRED + IPRED*EPS(1)

$THETA
(0, 992.275)   ; 1 Vm_pop
(0, 330.3)   ; 2 Km_pop
(0, 66.6863)   ; 3 V_pop

$OMEGA
0.100000148   ; ETA(1) Vm  (omega = 0.316228)
0.1400001956   ; ETA(2) Km  (omega = 0.374166)
0.100000148   ; ETA(3) V  (omega = 0.316228)

$SIGMA
0.04   ; EPS(1) DV proportional  (b = 0.2)

; -2LL: importance sampling (EONLY=1), proposal t with 5 df; first IMP adapts the proposal,
; the second evaluates with 10000 samples per individual (as Monolix). Reported value:
; 'OBJECTIVE FUNCTION VALUE WITH CONSTANT' of the last step.
$EST METHOD=COND INTERACTION FAST MAXEVAL=9999 NSIG=2 SIGL=6 PRINT=5 NOABORT
     FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=5 ISAMPLE=1000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=23456 NOCOV=1 FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=1 ISAMPLE=10000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=34567 NOCOV=1 FNLETA=0
