$PROBLEM Schoemaker oral1_1cpt_kaVVmKm | sch:oral1_1cpt_kaVVmKm runN030 | FOCEI
$ABBR DERIV2=NO
$INPUT ID TIME AMT DV MDV EVID CMT
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=tolswitch.f90
$MODEL COMP=(DEPOT,DEFDOSE) COMP=(CENTRAL,DEFOBS)

$PK
; individual parameters: theta_i = exp(MU + ETA), MU = log(typical value) [+ covariate effects];
; parameters without IIV have an ETA with OMEGA 0 FIX (MU referencing for FOCE FAST)
MU_1 = LOG(THETA(1))
VM = EXP(MU_1 + ETA(1))
MU_2 = LOG(THETA(2))
KM = EXP(MU_2 + ETA(2))
MU_3 = LOG(THETA(3))
V = EXP(MU_3 + ETA(3))
MU_4 = LOG(THETA(4))
KA = EXP(MU_4 + ETA(4))

$DES
CC = A(2)/V
DADT(1) = -KA*A(1)
DADT(2) = KA*A(1) - VM*CC/(KM + CC)

$ERROR
IPRED = A(2)/V
Y = IPRED + IPRED*EPS(1)

$THETA
(0, 1096.63)   ; 1 Vm_pop
(0, 298.867)   ; 2 Km_pop
(0, 90.0171)   ; 3 V_pop
(0, 1.2214)   ; 4 ka_pop

$OMEGA
0.1499997408   ; ETA(1) Vm  (omega = 0.387298)
0.1499997408   ; ETA(2) Km  (omega = 0.387298)
0.1499997408   ; ETA(3) V  (omega = 0.387298)
0.1499997408   ; ETA(4) ka  (omega = 0.387298)

$SIGMA
0.09   ; EPS(1) DV proportional  (b = 0.3)

; -2LL: importance sampling (EONLY=1), proposal t with 5 df; first IMP adapts the proposal,
; the second evaluates with 10000 samples per individual (as Monolix). Reported value:
; 'OBJECTIVE FUNCTION VALUE WITH CONSTANT' of the last step.
$EST METHOD=COND INTERACTION FAST MAXEVAL=9999 NSIG=2 SIGL=6 PRINT=5 NOABORT
     FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=5 ISAMPLE=1000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=23456 NOCOV=1 FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=1 ISAMPLE=10000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=34567 NOCOV=1 FNLETA=0
