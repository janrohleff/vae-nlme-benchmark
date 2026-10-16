$PROBLEM Myelosuppression | friberg rep03 | FOCEI
$ABBR DERIV2=NO
$INPUT ID TIME AMT DV MDV EVID CMT DVID
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=tolswitch.f90
$MODEL COMP=(CENTRAL,DEFDOSE,DEFOBS) COMP=(PROL) COMP=(TRANS1) COMP=(TRANS2) COMP=(TRANS3) COMP=(CIRC)

$PK
; individual parameters: theta_i = exp(MU + ETA), MU = log(typical value) [+ covariate effects];
; parameters without IIV have an ETA with OMEGA 0 FIX (MU referencing for FOCE FAST)
MU_1 = LOG(THETA(1))
V = EXP(MU_1 + ETA(1))
MU_2 = LOG(THETA(2))
CL = EXP(MU_2 + ETA(2))
MU_3 = LOG(THETA(3))
CIRC0 = EXP(MU_3 + ETA(3))
MU_4 = LOG(THETA(4))
MTT = EXP(MU_4 + ETA(4))
MU_5 = LOG(THETA(5))
GAM = EXP(MU_5 + ETA(5))
MU_6 = LOG(THETA(6))
SLOPE = EXP(MU_6 + ETA(6))
KTR = 4.0/MTT
A_0(2) = CIRC0
A_0(3) = CIRC0
A_0(4) = CIRC0
A_0(5) = CIRC0
A_0(6) = CIRC0

$DES
CC = A(1)/V
EDRUG = SLOPE*CC
DADT(1) = -CL/V*A(1)
DADT(2) = KTR*A(2)*(1 - EDRUG)*(CIRC0/A(6))**GAM - KTR*A(2)
DADT(3) = KTR*(A(2) - A(3))
DADT(4) = KTR*(A(3) - A(4))
DADT(5) = KTR*(A(4) - A(5))
DADT(6) = KTR*(A(5) - A(6))

$ERROR
IPRED = A(1)/V
IF (DVID.EQ.2) IPRED = A(6)
Y = IPRED + IPRED*EPS(1)
IF (DVID.EQ.2) Y = IPRED + IPRED*EPS(2)

$THETA
(0, 78.0)   ; 1 V_pop
(0, 21.0)   ; 2 Cl_pop
(0, 6.0)   ; 3 Circ0_pop
(0, 63.0)   ; 4 MTT_pop
(0, 0.255)   ; 5 gamma_pop
(0, 6.0)   ; 6 Slope_pop

$OMEGA
0.25   ; ETA(1) V  (omega = 0.5)
0.25   ; ETA(2) Cl  (omega = 0.5)
0.25   ; ETA(3) Circ0  (omega = 0.5)
0.25   ; ETA(4) MTT  (omega = 0.5)
0 FIX   ; ETA(5) gamma  (no IIV)
0.25   ; ETA(6) Slope  (omega = 0.5)

$SIGMA
0.09   ; EPS(1) Conc proportional  (b1 = 0.3)
0.09   ; EPS(2) ANC proportional  (b2 = 0.3)

; -2LL: importance sampling (EONLY=1), proposal t with 5 df; first IMP adapts the proposal,
; the second evaluates with 10000 samples per individual (as Monolix). Reported value:
; 'OBJECTIVE FUNCTION VALUE WITH CONSTANT' of the last step.
$EST METHOD=COND INTERACTION FAST MAXEVAL=9999 NSIG=2 SIGL=6 PRINT=5 NOABORT
     FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=5 ISAMPLE=1000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=23456 NOCOV=1 FORMAT=s1PE17.10
$EST METHOD=IMP INTERACTION EONLY=1 NITER=1 ISAMPLE=10000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=3 PRINT=1 SEED=34567 NOCOV=1 FNLETA=0
