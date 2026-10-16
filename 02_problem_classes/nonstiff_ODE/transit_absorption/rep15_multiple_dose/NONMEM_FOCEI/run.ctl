$PROBLEM Transit absorption | transit rep15 | FOCEI
$ABBR DERIV2=NO
$INPUT ID TIME AMT DV MDV EVID CMT
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=tolswitch.f90
$MODEL COMP=(TRANS1,DEFDOSE) COMP=(TRANS2) COMP=(TRANS3) COMP=(ABS) COMP=(CENTRAL,DEFOBS)

$PK
; individual parameters: theta_i = exp(MU + ETA), MU = log(typical value) [+ covariate effects];
; parameters without IIV have an ETA with OMEGA 0 FIX (MU referencing for FOCE FAST)
MU_1 = LOG(THETA(1))
MTT = EXP(MU_1 + ETA(1))
MU_2 = LOG(THETA(2))
KA = EXP(MU_2 + ETA(2))
MU_3 = LOG(THETA(3))
V = EXP(MU_3 + ETA(3))
MU_4 = LOG(THETA(4))
CL = EXP(MU_4 + ETA(4))
KTR = 4.0/MTT

$DES
DADT(1) = -KTR*A(1)
DADT(2) = KTR*(A(1) - A(2))
DADT(3) = KTR*(A(2) - A(3))
DADT(4) = KTR*A(3) - KA*A(4)
DADT(5) = KA*A(4) - CL/V*A(5)

$ERROR
IPRED = A(5)/V
Y = IPRED + IPRED*EPS(1)

$THETA
(0, 0.9)   ; 1 MTT_pop
(0, 1.05)   ; 2 ka_pop
(0, 70.0)   ; 3 V_pop
(0, 7.0)   ; 4 Cl_pop

$OMEGA
0.25   ; ETA(1) MTT  (omega = 0.5)
0.25   ; ETA(2) ka  (omega = 0.5)
0.25   ; ETA(3) V  (omega = 0.5)
0.25   ; ETA(4) Cl  (omega = 0.5)

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
