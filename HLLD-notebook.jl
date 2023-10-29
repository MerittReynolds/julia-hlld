### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 466b3e0a-f1ed-11ed-20b3-6122d6f2ba5e
"Given state vector U, returns primitive variables rho, vx, vy, vz, p, Bx, By, Bz"
function primitive_vars_from_state(U)
    
    # conservative variables
    rho = U[1]
    mx = U[2]
    my = U[3]
    mz = U[4]
    e = U[5]
    Bx = U[6]
    By = U[7]
    Bz = U[8]
    
    vx = mx/rho
    vy = my/rho
    vz = mz/rho

    gam = 5/3

    BB = Bx^2 + By^2 + Bz^2
    vv = vx^2 + vy^2 + vz^2
    p = (gam - 1)*(e - BB/2 - rho*vv/2)

    return rho, vx, vy, vz, p, Bx, By, Bz

end


# ╔═╡ 427ef6be-f15c-494c-975c-c9d9fe03f4aa
"Given primitive variables rho, vx, vy, vz, p, Bx, By, Bz returns state vector U"
function state_from_primitive_vars(rho, vx, vy, vz, p, Bx, By, Bz)
    
    mx = rho*vx
    my = rho*vy
    mz = rho*vz

    gam = 5/3

    BB = Bx^2 + By^2 + Bz^2
    vv = vx^2 + vy^2 + vz^2
    e = p/(gam - 1) + rho*vv/2 + BB/2

    return [rho, mx, my, mz, e, Bx, By, Bz]

end

# ╔═╡ 16de09cc-50d3-4ea3-b169-56c547ac4a96
"Given state vector U returns total energy e"
function totalenergy(U)
    
    return U[5]

end


# ╔═╡ 109374ec-707d-4e22-beeb-7575d981fbbf
"Given state vector U returns fast wave speed in x direction"
function cfast(U)
    
    #'fast wave speed in x direction'
    
    rho, vx, vy, vz, p, Bx, By, Bz = primitive_vars_from_state(U)
    
    gam = 5/3
    
    BB = Bx^2 + By^2 + Bz^2
    
    return sqrt((gam*p + BB + sqrt((gam*p + BB)^2 -4*gam*p*Bx^2))/(2*rho))

end


# ╔═╡ 65abed85-b95b-4254-9a83-250d86f11c9b
"Given state vector U returns total pressure (thermal plus magnetic)"
function ptotal(U)
    
    #'total pressure'
    
    rho, vx, vy, vz, p, Bx, By, Bz = primitive_vars_from_state(U)
    
    BB = Bx^2 + By^2 + Bz^2
    
    return p + BB/2

end

# ╔═╡ 68a08bae-406b-40fe-a456-6b2239ec32b4
"Calculate MHD fluxes"
function MHD_flux(U)
    
    #'MHD flux'

	e = totalenergy(U)
    rho, vx, vy, vz, p, Bx, By, Bz = primitive_vars_from_state(U)

    BB = Bx^2 + By^2 + Bz^2
    pT = p + BB/2

    F_rho = rho*vx
    F_mx = rho*vx^2 + pT - Bx^2
    F_my = rho*vy*vx - Bx*By
    F_mz = rho*vz*vx - Bx*Bz
    F_e = (e + pT)*vx - Bx*(vx*Bx + vy*By + vz*Bz)
    F_Bx = 0
    F_By = By*vx - Bx*vy
    F_Bz = Bz*vx - Bx*vz
    
    return [F_rho, F_mx, F_my, F_mz, F_e, F_Bx, F_By, F_Bz]

end

# ╔═╡ a079d247-7e87-4c99-8955-dbfd0823310a
"Calculate HLLD fluxes
    
Ref:
A multi-state HLL approximate Riemann solver for ideal magnetohydrodynamics
T. Miyoshi, K. Kusano
Journal of Computational Physics, 208 (2005) 315-344."
function HLLD_flux(U_L, U_R)
    
    # Calculate HLLD fluxes
    
    # Ref:
    # A multi-state HLL approximate Riemann solver for ideal
    # magnetohydrodynamics
    # T. Miyoshi, K. Kusano
    # Journal of Computational Physics, 208 (2005) 315-344.
    
    # Get necessary conservative variable total energy
    
    eL = totalenergy(U_L)
    eR = totalenergy(U_R)
    
    # Get left and right states as primitive variables
    
    rhoL, uL, vL, wL, pL, BxL, ByL, BzL = primitive_vars_from_state(U_L)
    rhoR, uR, vR, wR, pR, BxR, ByR, BzR = primitive_vars_from_state(U_R)
    
    # Normal magnetic field (should be continuous)
    # This is just in case the caller passes a discontinuous state
    
    Bx = (BxL + BxR)/2
    
    # Fluxes
    
    F_L = MHD_flux(U_L)
    F_R = MHD_flux(U_R)
    
    # Wave speeds
    
    cfL = cfast(U_L)
    cfR = cfast(U_R)
    
    # Total pressures
    
    pTL = ptotal(U_L)
    pTR = ptotal(U_R)
    
    # Signal speed estimation
    
    SL = min(uL, uR) - max(cfL, cfR)
    SR = max(uL, uR) + max(cfL, cfR)
   
    # Often-used velocity difference
    
    SmuR = SR - uR
    SmuL = SL - uL
    
    # Denominator for SM and pT expressions
    
    denSM = SmuR*rhoR - SmuL*rhoL
    
    # Average normal velocity in the Riemann fan
    # Equation (38)
    
    SM = (SmuR*rhoR*uR - SmuL*rhoL*uL - pTR + pTL)/denSM
    
    # Normal velocity is assumed to be constant over
    # the Riemann fan
    
    suL = SM
    ssuL = SM
    ssuR = SM
    suR = SM
    
    # Average total pressure in the Riemann fan
    # Miyoshi and Kusano equation (41)
    
    spT = (SmuR*rhoR*pTL - SmuL*rhoL*pTR
           + rhoL*rhoR*SmuR*SmuL*(uR - uL))/denSM
    spTL = spT
    spTR = spT
    
    # Density in the star regions
    
    srhoL = rhoL*SmuL/(SL - SM)
    srhoR = rhoR*SmuR/(SR - SM)
    
    # Denominator
    
    denL = rhoL*SmuL*(SL - SM) - Bx^2
    denR = rhoR*SmuR*(SR - SM) - Bx^2
    
    # Calculate v, By, w, Bz in star regions
    
    if abs(denL) > eps(1.0)
        svL = vL - Bx*ByL*(SM - uL)/denL
        sByL = ByL*(rhoL*SmuL^2 - Bx^2)/denL
        swL = wL - Bx*BzL*(SM - uL)/denL
        sBzL = BzL*(rhoL*SmuL^2 - Bx^2)/denL
    else
        svL = vL
        sByL = ByL
        swL = wL
        sBzL = BzL
	end
        
    if abs(denR) > eps(1.0)
        svR = vR - Bx*ByR*(SM - uR)/denR
        sByR = ByR*(rhoR*SmuR^2 - Bx^2)/denR
        swR = wR - Bx*BzR*(SM - uR)/denR
        sBzR = BzR*(rhoR*SmuR^2 - Bx^2)/denR
    else
        svR = vR
        sByR = ByR
        swR = wR
        sBzR = BzR
	end
    
    # Calculate energy density in the star regions
    # Miyoshi and Kusano equation (48)

    vBL = uL*Bx + vL*ByL + wL*BzL
    svBL = suL*Bx + svL*sByL + swL*sBzL
    seL = (SmuL*eL - pTL*uL + spT*SM + Bx*(vBL - svBL))/(SL - SM)

    vBR = uR*Bx + vR*ByR + wR*BzR
    svBR = suR*Bx + svR*sByR + swR*sBzR
    seR = (SmuR*eR - pTR*uR + spT*SM + Bx*(vBR - svBR))/(SR - SM)
    
    # Propagation speeds of Alfven wave in intermediate states
    # Equation (51)
    
    sqrt_srhoL = sqrt(srhoL)
    sqrt_srhoR = sqrt(srhoR)    
    
    sSL = SM - abs(Bx)/sqrt_srhoL
    sSR = SM + abs(Bx)/sqrt_srhoR

    # Density and pressure in the star-star regions
    
    ssrhoL = srhoL
    ssrhoR = srhoR
    
    sspTL = spTL
    sspTR = spTR

    # Equations (59)-(62)

    sign_Bx = sign(Bx)
    
    sqrt_srhoLsrhoR = sqrt_srhoL*sqrt_srhoR
    
    den = sqrt_srhoL + sqrt_srhoR
    
    ssv = (sqrt_srhoL*svL + sqrt_srhoR*svR + (sByR - sByL)*sign_Bx)/den
    ssw = (sqrt_srhoL*swL + sqrt_srhoR*swR + (sBzR - sBzL)*sign_Bx)/den
    ssBy = (sqrt_srhoL*sByR + sqrt_srhoR*sByL 
            + sqrt_srhoLsrhoR*(svR - svL)*sign_Bx)/den
    ssBz = (sqrt_srhoL*sBzR + sqrt_srhoR*sBzL 
            + sqrt_srhoLsrhoR*(swR - swL)*sign_Bx)/den
    
    # TODO: double check these assumptions
    
    ssvL = ssv
    ssvR = ssv
    
    sswL = ssw
    sswR = ssw
    
    ssByL = ssBy
    ssByR = ssBy
    
    ssBzL = ssBz
    ssBzR = ssBz
    
    # Energy density in star-star regions
    
    ssvBL = ssuL*Bx + ssvL*ssByL + sswL*ssBzL
    ssvBR = ssuR*Bx + ssvR*ssByR + sswR*ssBzR
    
    sseL = seL - sqrt_srhoL*(svBL - ssvBL)*sign_Bx
    sseR = seR + sqrt_srhoR*(svBR - ssvBR)*sign_Bx

    # States (conservative variables) as arrays so we can use in flux expressions
    
    sU_L = [srhoL, srhoL*suL, srhoL*svL, srhoL*swL, seL, Bx, sByL, sBzL]
    sU_R = [srhoR, srhoR*suR, srhoR*svR, srhoR*swR, seR, Bx, sByR, sBzR]

    ssU_L = [ssrhoL, ssrhoL*ssuL, ssrhoL*ssvL, ssrhoL*sswL, sseL, Bx, ssByL, ssBzL]
    ssU_R = [ssrhoR, ssrhoR*ssuR, ssrhoR*ssvR, ssrhoR*sswR, sseR, Bx, ssByR, ssBzR]

    # Fluxes
    
    sF_L = F_L + SL*sU_L - SL*U_L
    ssF_L = F_L + sSL*ssU_L - (sSL - SL)*sU_L - SL*U_L

    sF_R = F_R + SR*sU_R - SR*U_R
    ssF_R = F_R + sSR*ssU_R - (sSR - SR)*sU_R - SL*U_R
        
    if SL > 0
        F_HLLD = F_L
	elseif sSL >= 0
        F_HLLD = sF_L
    elseif SM >= 0
        F_HLLD = ssF_L
    elseif sSR >= 0
        F_HLLD = ssF_R
    elseif SR >= 0
        F_HLLD = sF_R
    else
        F_HLLD = F_R
	end
    
    return F_HLLD
end

# ╔═╡ 438a6c65-2ccd-428d-82c1-b3a52778652e
begin
	# Compare HLLD and MHD fluxes
	rho = 1
	vx = 1
	vy = 0.5
	vz = 0
	Bx = 0.001
	By = 0.3
	Bz = 0.1
	p = 1

	U = state_from_primitive_vars(rho, vx, vy, vz, p, Bx, By, Bz)

	println("Compare HLLD and MHD fluxes")
	println("   U: ", U)

	F_MHD = MHD_flux(U)
	F = HLLD_flux(U, U)

	println("HLLD: ", F)
	println(" MHD: ", F_MHD)

end


# ╔═╡ c36f10c7-98a9-4ffb-b2c2-ce97ce257e91
begin
# See if we can sustain a current sheet

rhoL = 1
vxL = 0
vyL = 0
vzL = 0
BxL = 0
ByL = 0
BzL = 0
pL = 1.5

U_L = state_from_primitive_vars(rhoL, vxL, vyL, vzL, pL, BxL, ByL, BzL)

rhoR = 1
vxR = 0
vyR = 0
vzR = 0
BxR = 0
ByR = 1
BzR = 0
pR = 1

U_R = state_from_primitive_vars(rhoR, vxR, vyR, vzR, pR, BxR, ByR, BzR)

# HLLD flux
xF = HLLD_flux(U_L, U_R)

	println("See if we can sustain a current sheet")
	println("HLLD: ", xF)

end

# ╔═╡ bc77efb7-42d1-4d22-b556-172546cc4c8d


# ╔═╡ Cell order:
# ╟─466b3e0a-f1ed-11ed-20b3-6122d6f2ba5e
# ╟─427ef6be-f15c-494c-975c-c9d9fe03f4aa
# ╟─16de09cc-50d3-4ea3-b169-56c547ac4a96
# ╟─109374ec-707d-4e22-beeb-7575d981fbbf
# ╟─65abed85-b95b-4254-9a83-250d86f11c9b
# ╠═68a08bae-406b-40fe-a456-6b2239ec32b4
# ╟─a079d247-7e87-4c99-8955-dbfd0823310a
# ╟─438a6c65-2ccd-428d-82c1-b3a52778652e
# ╟─c36f10c7-98a9-4ffb-b2c2-ce97ce257e91
# ╠═bc77efb7-42d1-4d22-b556-172546cc4c8d
