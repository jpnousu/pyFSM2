#-----------------------------------------------------------------------
# Properties of vegetation canopy layers
#-----------------------------------------------------------------------

use CONSTANTS, only: &
  hcap_ice            ! Specific heat capacity of ice (J/K/kg)

use LAYERS, only: &
  Ncnpy,             &! Number of canopy layers
  fvg1                ! Fraction of vegetation in upper canopy layer

use PARAMETERS, only: &
  cvai,              &! Vegetation heat capacity per unit VAI (J/K/m^2)
  svai                ! Intercepted snow capacity per unit VAI (kg/m^2)

implicit none

real, intent(in) :: &
  Sveg(Ncnpy),       &! Snow mass on vegetation layers (kg/m^2)
  Tveg(Ncnpy),       &! Vegetation layer temperatures (K)
  VAI                 ! Vegetation area index

real, intent(out) :: &
  cveg(Ncnpy),       &! Vegetation heat capacities (J/K/m^2)
  fcans(Ncnpy),      &! Canopy layer snowcover fractions
  lveg(Ncnpy),       &! Canopy layer vegetation area indices
  Scap(Ncnpy),       &! Canopy layer snow capacities (kg/m^2)
  Tveg0(Ncnpy)        ! Vegetation temperatures at start of timestep (K)

#if CANMOD == 1
lveg(1) = VAI
#endif
#if CANMOD == 2
lveg(1) = fvg1*VAI
lveg(2) = (1 - fvg1)*VAI
#endif
cveg(:) = cvai*lveg(:) + hcap_ice*Sveg(:)
Scap = svai*lveg(:)
fcans(:) = 0
if (VAI > epsilon(VAI)) fcans(:) = (Sveg(:)/Scap(:))**0.67
Tveg0 = Tveg

end subroutine CANOPY