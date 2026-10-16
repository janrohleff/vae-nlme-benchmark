$PROBLEM Tumour growth inhibition | tgi rep18 | FOCEI
$ABBR DERIV2=NO
$INPUT ID TIME AMT DV MDV EVID CMT
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=3 ATOL=6
$MODEL COMP=(PKCMT,DEFDOSE) COMP=(X1,DEFOBS) COMP=(X2) COMP=(X3) COMP=(X4)

$PK
; individual parameters: theta_i = exp(MU + ETA), MU = log(typical value) [+ covariate effects];
; parameters without IIV have an ETA with OMEGA 0 FIX (MU referencing for FOCE FAST)
KEL = 5.0   ; fixed, as in Monolix and VAE Studio
V = 1.0   ; fixed, as in Monolix and VAE Studio
MU_1 = LOG(THETA(1))
LAM0 = EXP(MU_1 + ETA(1))
MU_2 = LOG(THETA(2))
LAM1 = EXP(MU_2 + ETA(2))
MU_3 = LOG(THETA(3))
K1 = EXP(MU_3 + ETA(3))
MU_4 = LOG(THETA(4))
K2 = EXP(MU_4 + ETA(4))
MU_5 = LOG(THETA(5))
W0 = EXP(MU_5 + ETA(5))
A_0(2) = W0

$DES
CC = A(1)/V
WT = A(2) + A(3) + A(4) + A(5)
DADT(1) = -KEL*A(1)
DADT(2) = LAM0*A(2)/(1 + (LAM0/LAM1*WT)**20)**0.05 - K2*CC*A(2)
DADT(3) = K2*CC*A(2) - K1*A(3)
DADT(4) = K1*(A(3) - A(4))
DADT(5) = K1*(A(4) - A(5))

$ERROR
IPRED = A(2) + A(3) + A(4) + A(5)
Y = IPRED + IPRED*EPS(1)

$THETA
(0, 0.175)   ; 1 lambda0_pop
(0, 0.225)   ; 2 lambda1_pop
(0, 0.6)   ; 3 k1_pop
(0, 0.75)   ; 4 k2_pop
(0, 0.006)   ; 5 w0_pop

$OMEGA
0.25   ; ETA(1) lambda0  (omega = 0.5)
0.25   ; ETA(2) lambda1  (omega = 0.5)
0 FIX   ; ETA(3) k1  (no IIV)
0.25   ; ETA(4) k2  (omega = 0.5)
0.25   ; ETA(5) w0  (omega = 0.5)

$SIGMA
0.09   ; EPS(1) Tumor proportional  (b = 0.3)

; -2LL: importance sampling (EONLY=1), proposal t with 5 df; first IMP adapts the proposal,
; the second evaluates with 10000 samples per individual (as Monolix). Reported value:
; 'OBJECTIVE FUNCTION VALUE WITH CONSTANT' of the last step.
$EST METHOD=COND INTERACTION FAST MAXEVAL=9999 NSIG=1 SIGL=3 PRINT=5 NOABORT
     FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=5 ISAMPLE=1000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=23456 NOCOV=1 FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=1 ISAMPLE=10000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=34567 NOCOV=1 FNLETA=0
