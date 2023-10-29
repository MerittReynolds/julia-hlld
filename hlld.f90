! Calculate HLLD fluxes
!
! Ref:
! A multi-state HLL approximate Riemann solver for ideal magnetohydrodynamics
! T. Miyoshi, K. Kusano
! Journal of Computational Physics, 208 (2005) 315-344."

module hlld

integer,parameter :: dp = kind(0.d0)

real(dp),parameter :: gam = 5.0_dp/3.0_dp

private dp
private gam

contains

!==============================================================================
pure subroutine primitive_vars_from_state(U, rho, vx, vy, vz, p, Bx, By, Bz)
    
implicit real(dp) (a-z)

real(dp),dimension(8),intent(in) :: U
real(dp),intent(out) :: rho, vx, vy, vz, p, Bx, By, Bz
!------------------------------------------------------------------------------

! conservative variables
rho = U(1)
mx = U(2)
my = U(3)
mz = U(4)
e = U(5)
Bx = U(6)
By = U(7)
Bz = U(8)
    
! velocity from momentum
vx = mx/rho
vy = my/rho
vz = mz/rho

! pressure from total energy
BB = Bx**2 + By**2 + Bz**2
vv = vx**2 + vy**2 + vz**2
p = (gam - 1)*(e - BB/2 - rho*vv/2)

end subroutine primitive_vars_from_state

!==============================================================================
pure subroutine state_from_primitive_vars(rho, vx, vy, vz, p, Bx, By, Bz, U)
    
implicit real(dp) (a-z)

real(dp),intent(in) :: rho, vx, vy, vz, p, Bx, By, Bz
real(dp),dimension(8),intent(out) :: U
!------------------------------------------------------------------------------

mx = rho*vx
my = rho*vy
mz = rho*vz

BB = Bx**2 + By**2 + Bz**2
vv = vx**2 + vy**2 + vz**2
e = p/(gam - 1) + rho*vv/2 + BB/2

U(1) = rho
U(2) = mx
U(3) = my
U(4) = mz
U(5) = e
U(6) = Bx
U(7) = By
U(8) = Bz

end subroutine state_from_primitive_vars

!==============================================================================
pure function totalenergy(U)

real(dp),dimension(8),intent(in) :: U
real(dp) :: totalenergy
!------------------------------------------------------------------------------

totalenergy = U(5)

end function totalenergy

!==============================================================================
pure function cfast(U)
!'fast wave speed in x direction'

implicit real(dp) (a-z)
real(dp),dimension(8),intent(in) :: U
real(dp) :: cfast
!------------------------------------------------------------------------------

call primitive_vars_from_state(U, rho, vx, vy, vz, p, Bx, By, Bz)
    
BB = Bx**2 + By**2 + Bz**2
    
cfast = sqrt((gam*p + BB + sqrt((gam*p + BB)**2 -4*gam*p*Bx**2))/(2*rho))

end function cfast

!==============================================================================
!'total pressure'
pure function ptotal(U)
    
implicit real(dp) (a-z)
real(dp),dimension(8),intent(in) :: U
real(dp) :: ptotal
!------------------------------------------------------------------------------
    
call primitive_vars_from_state(U, rho, vx, vy, vz, p, Bx, By, Bz)
    
BB = Bx**2 + By**2 + Bz**2
    
ptotal = p + BB/2

end function ptotal

!==============================================================================
pure function MHD_flux(U) result(F)

! MHD flux, returns array

implicit real(dp) (a-z)
real(dp),dimension(8),intent(in) :: U
real(dp),dimension(8) :: F
!------------------------------------------------------------------------------

e = totalenergy(U)
call primitive_vars_from_state(U, rho, vx, vy, vz, p, Bx, By, Bz)

BB = Bx**2 + By**2 + Bz**2
pT = p + BB/2

F_rho = rho*vx
F_mx = rho*vx**2 + pT - Bx**2
F_my = rho*vy*vx - Bx*By
F_mz = rho*vz*vx - Bx*Bz
F_e = (e + pT)*vx - Bx*(vx*Bx + vy*By + vz*Bz)
F_Bx = 0
F_By = By*vx - Bx*vy
F_Bz = Bz*vx - Bx*vz
    
F = [F_rho, F_mx, F_my, F_mz, F_e, F_Bx, F_By, F_Bz]

end function MHD_flux

!==============================================================================
pure function HLLD_flux(U_L, U_R) result(F_HLLD)

! Calculate HLLD fluxes
!
! Ref:
! A multi-state HLL approximate Riemann solver for ideal
! magnetohydrodynamics
! T. Miyoshi, K. Kusano
! Journal of Computational Physics, 208 (2005) 315-344.

implicit real(dp) (a-z)

real(dp),dimension(8),intent(in) :: U_L, U_R
real(dp),dimension(8) :: F_HLLD

real(dp),parameter :: one = 1.0_dp

real(dp),dimension(8) :: F_L, F_R, sF_L, sF_R, ssF_L, ssF_R
real(dp),dimension(8) :: sU_L, sU_r, ssU_L, ssU_R
!------------------------------------------------------------------------------


! Get necessary conservative variable total energy

eL = totalenergy(U_L)
eR = totalenergy(U_R)

! Get left and right states as primitive variables

call primitive_vars_from_state(U_L, rhoL, uL, vL, wL, pL, BxL, ByL, BzL)
call primitive_vars_from_state(U_R, rhoR, uR, vR, wR, pR, BxR, ByR, BzR)

! Normal magnetic field (should be continuous)
! This is just in case the caller passes a discontinuous state

Bx = (BxL + BxR)/2

! Fluxes

F_L = MHD_flux(U_L)
F_R = MHD_flux(U_R)

! Wave speeds

cfL = cfast(U_L)
cfR = cfast(U_R)

! Total pressures

pTL = ptotal(U_L)
pTR = ptotal(U_R)

! Signal speed estimation

SL = min(uL, uR) - max(cfL, cfR)
SR = max(uL, uR) + max(cfL, cfR)

! Often-used velocity difference

SmuR = SR - uR
SmuL = SL - uL

! Denominator for SM and pT expressions

denSM = SmuR*rhoR - SmuL*rhoL

! Average normal velocity in the Riemann fan
! Equation (38)

SM = (SmuR*rhoR*uR - SmuL*rhoL*uL - pTR + pTL)/denSM

! Normal velocity is assumed to be constant over
! the Riemann fan

suL = SM
ssuL = SM
ssuR = SM
suR = SM

! Average total pressure in the Riemann fan
! Miyoshi and Kusano equation (41)

spT = (SmuR*rhoR*pTL - SmuL*rhoL*pTR &
       + rhoL*rhoR*SmuR*SmuL*(uR - uL))/denSM
spTL = spT
spTR = spT

! Density in the star regions

srhoL = rhoL*SmuL/(SL - SM)
srhoR = rhoR*SmuR/(SR - SM)

! Denominator

denL = rhoL*SmuL*(SL - SM) - Bx**2
denR = rhoR*SmuR*(SR - SM) - Bx**2

! Calculate v, By, w, Bz in star regions

if(abs(denL) > epsilon(1.0))then
    svL = vL - Bx*ByL*(SM - uL)/denL
    sByL = ByL*(rhoL*SmuL**2 - Bx**2)/denL
    swL = wL - Bx*BzL*(SM - uL)/denL
    sBzL = BzL*(rhoL*SmuL**2 - Bx**2)/denL
else
    svL = vL
    sByL = ByL
    swL = wL
    sBzL = BzL
endif

if(abs(denR) > epsilon(1.0))then
    svR = vR - Bx*ByR*(SM - uR)/denR
    sByR = ByR*(rhoR*SmuR**2 - Bx**2)/denR
    swR = wR - Bx*BzR*(SM - uR)/denR
    sBzR = BzR*(rhoR*SmuR**2 - Bx**2)/denR
else
    svR = vR
    sByR = ByR
    swR = wR
    sBzR = BzR
endif

! Calculate energy density in the star regions
! Miyoshi and Kusano equation (48)

vBL = uL*Bx + vL*ByL + wL*BzL
svBL = suL*Bx + svL*sByL + swL*sBzL
seL = (SmuL*eL - pTL*uL + spT*SM + Bx*(vBL - svBL))/(SL - SM)

vBR = uR*Bx + vR*ByR + wR*BzR
svBR = suR*Bx + svR*sByR + swR*sBzR
seR = (SmuR*eR - pTR*uR + spT*SM + Bx*(vBR - svBR))/(SR - SM)

! Propagation speeds of Alfven wave in intermediate states
! Equation (51)

sqrt_srhoL = sqrt(srhoL)
sqrt_srhoR = sqrt(srhoR)

sSL = SM - abs(Bx)/sqrt_srhoL
sSR = SM + abs(Bx)/sqrt_srhoR

! Density and pressure in the star-star regions

ssrhoL = srhoL
ssrhoR = srhoR

sspTL = spTL
sspTR = spTR

! Equations (59)-(62)

sign_Bx = sign(one, Bx)

sqrt_srhoLsrhoR = sqrt_srhoL*sqrt_srhoR

den = sqrt_srhoL + sqrt_srhoR

ssv = (sqrt_srhoL*svL + sqrt_srhoR*svR + (sByR - sByL)*sign_Bx)/den
ssw = (sqrt_srhoL*swL + sqrt_srhoR*swR + (sBzR - sBzL)*sign_Bx)/den
ssBy = (sqrt_srhoL*sByR + sqrt_srhoR*sByL &
        + sqrt_srhoLsrhoR*(svR - svL)*sign_Bx)/den
ssBz = (sqrt_srhoL*sBzR + sqrt_srhoR*sBzL &
        + sqrt_srhoLsrhoR*(swR - swL)*sign_Bx)/den

! TODO: double check these assumptions

ssvL = ssv
ssvR = ssv

sswL = ssw
sswR = ssw

ssByL = ssBy
ssByR = ssBy

ssBzL = ssBz
ssBzR = ssBz

! Energy density in star-star regions

ssvBL = ssuL*Bx + ssvL*ssByL + sswL*ssBzL
ssvBR = ssuR*Bx + ssvR*ssByR + sswR*ssBzR

sseL = seL - sqrt_srhoL*(svBL - ssvBL)*sign_Bx
sseR = seR + sqrt_srhoR*(svBR - ssvBR)*sign_Bx

! States (conservative variables) as arrays so we can use in flux expressions

sU_L = [srhoL, srhoL*suL, srhoL*svL, srhoL*swL, seL, Bx, sByL, sBzL]
sU_R = [srhoR, srhoR*suR, srhoR*svR, srhoR*swR, seR, Bx, sByR, sBzR]

ssU_L = [ssrhoL, ssrhoL*ssuL, ssrhoL*ssvL, ssrhoL*sswL, sseL, Bx, ssByL, ssBzL]
ssU_R = [ssrhoR, ssrhoR*ssuR, ssrhoR*ssvR, ssrhoR*sswR, sseR, Bx, ssByR, ssBzR]

! Fluxes

sF_L = F_L + SL*sU_L - SL*U_L
ssF_L = F_L + sSL*ssU_L - (sSL - SL)*sU_L - SL*U_L

sF_R = F_R + SR*sU_R - SR*U_R
ssF_R = F_R + sSR*ssU_R - (sSR - SR)*sU_R - SL*U_R

if (SL > 0) then
    F_HLLD = F_L
elseif (sSL >= 0) then
    F_HLLD = sF_L
elseif (SM >= 0) then
    F_HLLD = ssF_L
elseif (sSR >= 0) then
    F_HLLD = ssF_R
elseif (SR >= 0) then
    F_HLLD = sF_R
else
    F_HLLD = F_R
endif

end function HLLD_flux

end module hlld
