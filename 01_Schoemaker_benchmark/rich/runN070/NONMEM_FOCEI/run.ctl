$PROBLEM Schoemaker oral1_2cpt_kaV1QV2VmKm | sch:oral1_2cpt_kaV1QV2VmKm runN070 | FOCEI
$ABBR DERIV2=NO
$INPUT ID TIME AMT DV MDV EVID CMT
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=tolswitch.f90
$MODEL COMP=(DEPOT,DEFDOSE) COMP=(CENTRAL,DEFOBS) COMP=(PERIPH)

$PK
; individual parameters: theta_i = exp(MU + ETA), MU = log(typical value) [+ covariate effects];
; parameters without IIV have an ETA with OMEGA 0 FIX (MU referencing for FOCE FAST)
MU_1 = LOG(THETA(1))
VM = EXP(MU_1 + ETA(1))
MU_2 = LOG(THETA(2))
KM = EXP(MU_2 + ETA(2))
MU_3 = LOG(THETA(3))
V1 = EXP(MU_3 + ETA(3))
MU_4 = LOG(THETA(4))
Q = EXP(MU_4 + ETA(4))
MU_5 = LOG(THETA(5))
V2 = EXP(MU_5 + ETA(5))
MU_6 = LOG(THETA(6))
KA = EXP(MU_6 + ETA(6))

$DES
CC = A(2)/V1
DADT(1) = -KA*A(1)
DADT(2) = KA*A(1) - VM*CC/(KM + CC) - Q/V1*A(2) + Q/V2*A(3)
DADT(3) = Q/V1*A(2) - Q/V2*A(3)

$ERROR
IPRED = A(2)/V1
Y = IPRED + IPRED*EPS(1)

$THETA
(0, 1211.97)   ; 1 Vm_pop
(0, 298.867)   ; 2 Km_pop
(0, 90.0171)   ; 3 V1_pop
(0, 4.95303)   ; 4 Q_pop
(0, 60.3403)   ; 5 V2_pop
(0, 1.24608)   ; 6 ka_pop

$OMEGA
0.1499997408   ; ETA(1) Vm  (omega = 0.387298)
0.1499997408   ; ETA(2) Km  (omega = 0.387298)
0.1499997408   ; ETA(3) V1  (omega = 0.387298)
0.1499997408   ; ETA(4) Q  (omega = 0.387298)
0.1499997408   ; ETA(5) V2  (omega = 0.387298)
0.1499997408   ; ETA(6) ka  (omega = 0.387298)

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
