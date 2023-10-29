program main

use hlld

implicit none

integer,parameter :: dp = kind(0.d0)

! For test 1
real(dp) :: rho, vx, vy, vz, Bx, By, Bz, p
real(dp),dimension(8) :: U, F, F_MHD

! For test 2
real(dp) :: rhoL, vxL, vyL, vzL, BxL, ByL, BzL, pL
real(dp) :: rhoR, vxR, vyR, vzR, BxR, ByR, BzR, pR
real(dp),dimension(8) :: U_L, U_R, xF

! Compare MHD and HLLD flux

print *, 'TEST 1 - comparison of MHD and HLLD fluxes'

rho = 1
vx = 1
vy = 0.5
vz = 0
Bx = 0.001
By = 0.3
Bz = 0.1
p = 1

call state_from_primitive_vars(rho, vx, vy, vz, p, Bx, By, Bz, U)

write(*,'(A,8(F12.5,1X))') 'U:    ', U

F_MHD = MHD_flux(U)
F = HLLD_flux(U, U)

write(*,'(A,8(F12.5,1X))') 'HLLD: ', F
write(*,'(A,8(F12.5,1X))') 'MHD:  ', F_MHD
write(*,'(A,8(F12.5,1X))') 'diff: ', F - F_MHD

! See if we can sustain a current sheet

print *
print *, 'TEST 2 - current sheet'

rhoL = 1
vxL = 0
vyL = 0
vzL = 0
BxL = 0
ByL = 0
BzL = 0
pL = 1.5

call state_from_primitive_vars(rhoL, vxL, vyL, vzL, pL, BxL, ByL, BzL, U_L)

rhoR = 1
vxR = 0
vyR = 0
vzR = 0
BxR = 0
ByR = 1
BzR = 0
pR = 1

call state_from_primitive_vars(rhoR, vxR, vyR, vzR, pR, BxR, ByR, BzR, U_R)

! HLLD flux
xF = HLLD_flux(U_L, U_R)

write(*,'(A,8(F12.5,1X))') 'U_L:  ', U_L
write(*,'(A,8(F12.5,1X))') 'U_R:  ', U_R
write(*,'(A,8(F12.5,1X))') 'HLLD: ', xF

end program main
