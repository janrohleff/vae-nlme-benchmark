$PROBLEM Schoemaker infusion_1cpt_VCl | sch:infusion_1cpt_VCl runN012 | FOCEI
$ABBR DERIV2=NO
$INPUT ID TIME AMT RATE DV MDV EVID CMT
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN1 TRANS2

$PK
; individual parameters: theta_i = exp(MU + ETA), MU = log(typical value) [+ covariate effects];
; parameters without IIV have an ETA with OMEGA 0 FIX (MU referencing for FOCE FAST)
MU_1 = LOG(THETA(1))
CL = EXP(MU_1 + ETA(1))
MU_2 = LOG(THETA(2))
V = EXP(MU_2 + ETA(2))
S1 = V

$ERROR
IPRED = F
Y = IPRED + IPRED*EPS(1)

$THETA
(0, 4.95303)   ; 1 Cl_pop
(0, 90.0171)   ; 2 V_pop

$OMEGA
0.1499997408   ; ETA(1) Cl  (omega = 0.387298)
0.1499997408   ; ETA(2) V  (omega = 0.387298)

$SIGMA
0.09   ; EPS(1) DV proportional  (b = 0.3)

; -2LL: importance sampling (EONLY=1), proposal t with 5 df; first IMP adapts the proposal,
; the second evaluates with 10000 samples per individual (as Monolix). Reported value:
; 'OBJECTIVE FUNCTION VALUE WITH CONSTANT' of the last step.
$EST METHOD=COND INTERACTION FAST MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT
     FORMAT=s1PE17.10 PARAFILE=OFF
$EST METHOD=IMP INTERACTION EONLY=1 NITER=5 ISAMPLE=1000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=9 PRINT=1 SEED=23456 NOCOV=1 FORMAT=s1PE17.10 PARAFILE=ON
$EST METHOD=IMP INTERACTION EONLY=1 NITER=1 ISAMPLE=10000 DF=5 ISCALE_MIN=1 ISCALE_MAX=1 MAPITER=0
     CTYPE=0 SIGL=9 PRINT=1 SEED=34567 NOCOV=1 FNLETA=0 PARAFILE=ON
