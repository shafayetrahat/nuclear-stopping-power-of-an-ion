# nuclear-stopping-power-of-an-ion

This project investigates the nuclear stopping power of ions penetrating matter using two models: the
Ziegler, Biersack, Littmark (ZBL) screened Coulomb potential and a universal empirical formula. The
stopping power, defined as the energy loss per unit path length of an ion, was computed numerically for
two ion-target combinations: hydrogen (H) in silicon (Si) and silicon (Si) in gold (Au). The implemen-
tation relied on numerical integration (Gauss-Legendre quadrature) and root-finding (bisection method)
to solve the scattering integrals. Results were compared between the ZBL and universal models, proving
good agreement at intermediate energies but deviations at low and high energies due to approximations
in the universal formula. This study demonstrates the effectiveness of numerical methods in simulating
ion-matter interactions while highlighting limitations in empirical models.

Details can be found in the report.pdf

# Relative difference in nuclear Stopping power of H --> Si and Si --> Au

![Nuclear_Stopping_Power_Comparison](https://github.com/user-attachments/assets/706ef16c-4476-471f-be8d-16d22140ae7d)

# Relative_Error_Nuclear_Stopping_Power

![Relative_Error_Nuclear_Stopping_Power](https://github.com/user-attachments/assets/72c69ff0-543e-445a-b0f0-b4288c2bd466)
