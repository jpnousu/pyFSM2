#-----------------------------------------------------------------------
# Thermal properties of snow and soil
#-----------------------------------------------------------------------

subroutine THERMAL(Dsnw,Nsnow,Sice,Sliq,Tsnow,Tsoil,Vsmc,              &
                   csoil,Ds1,gs1,ksnow,ksoil,ks1,Ts1)

! Thermal conductivity of snow
ksnow = kfix
#if CONDCT == 1
do k = 1, Nsnow
  rhos = rhof
#if DENSTY == 1
  if (Dsnw(k) > epsilon(Dsnw)) rhos = (Sice(k) + Sliq(k)) / Dsnw(k)
#endif
  ksnow(k) = 2.224*(rhos/rho_wat)**1.885
end do
#endif

! Heat capacity and thermal conductivity of soil
dPsidT = - rho_ice*Lf/(rho_wat*g*Tm)
do k = 1, Nsoil
  csoil(k) = hcap_soil*Dzsoil(k)
  ksoil(k) = hcon_soil
  if (Vsmc(k) > epsilon(Vsmc)) then
    dthudT = 0
    sthu = Vsmc(k)
    sthf = 0
    Tc = Tsoil(k) - Tm
    Tmax = Tm + (sathh/dPsidT)*(Vsat/Vsmc(k))**b
    if (Tsoil(k) < Tmax) then
      dthudT = (-dPsidT*Vsat/(b*sathh)) * (dPsidT*Tc/sathh)**(-1/b - 1)
      sthu = Vsat*(dPsidT*Tc/sathh)**(-1/b)
      sthu = min(sthu, Vsmc(k))
      sthf = (Vsmc(k) - sthu)*rho_wat/rho_ice
    end if
    Mf = rho_ice*Dzsoil(k)*sthf
    Mu = rho_wat*Dzsoil(k)*sthu
    csoil(k) = hcap_soil*Dzsoil(k) + hcap_ice*Mf + hcap_wat*Mu +       &
               rho_wat*Dzsoil(k)*((hcap_wat - hcap_ice)*Tc + Lf)*dthudT
    Smf = rho_ice*sthf/(rho_wat*Vsat)
    Smu = sthu/Vsat
    thice = 0
    if (Smf > 0) thice = Vsat*Smf/(Smu + Smf) 
    thwat = 0
    if (Smu > 0) thwat = Vsat*Smu/(Smu + Smf)
    hcon_sat = hcon_soil*(hcon_wat**thwat)*(hcon_ice**thice) /         &
              (hcon_air**Vsat)
    ksoil(k) = (hcon_sat - hcon_soil)*(Smf + Smu) + hcon_soil
    if (k == 1) gs1 = gsat*max((Smu*Vsat/Vcrit)**2, 1.)
  end if
end do

! Surface layer
Ds1 = max(Dzsoil(1), Dsnw(1))
Ts1 = Tsoil(1) + (Tsnow(1) - Tsoil(1))*Dsnw(1)/Dzsoil(1)
ks1 = Dzsoil(1)/(2*Dsnw(1)/ksnow(1) + (Dzsoil(1) - 2*Dsnw(1))/ksoil(1))
snd = sum(Dsnw)
if (snd > 0.5*Dzsoil(1)) ks1 = ksnow(1)
if (snd > Dzsoil(1)) Ts1 = Tsnow(1)

end subroutine THERMAL