# This script tests the convergence of the nuclear stopping power calculation
# by varying the quadrature order and b_max.
# It uses the Ziegler-Biersack-Littmark (ZBL) model for nuclear stopping power.
# The script calculates the stopping power for different projectile-target combinations
# and energies, and saves the results to a CSV file.
# The script uses the Gauss-Legendre quadrature method for numerical integration.
# Author: Md Shafayet Islam
# Date: 2025-05-22

import nuclear_stopping as nstp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def nuclear_stopping_ZBL_convergence_test(E_lab, Z1, Z2, M1, M2, b_max_factor, order_list):
    """
    Test convergence of nuclear stopping power calculation by varying quadrature order and b_max.
    Args:
        E_lab (float): Lab energy in eV.
        Z1, Z2 (int): Atomic numbers of projectile and target.
        M1, M2 (float): Masses of projectile and target (amu).
        b_max_factor (float): Multiplier for screening length (a_u) to set b_max.
        order_list (list): List of quadrature orders to test.
    Returns:
        dict: Results for different quadrature orders.
    """
    E_com = (M2 / (M1 + M2)) * E_lab
    a_u = 0.46848 / (Z1**0.23 + Z2**0.23)
    gamma = 4 * M1 * M2 / (M1 + M2)**2
    b_max = b_max_factor * a_u

    results = {}
    for order in order_list:
        # Gauss-Legendre quadrature for b integral
        nodes, weights = np.polynomial.legendre.leggauss(order)
        b_scaled = 0.5 * b_max * (nodes + 1)  # Transform to [0, b_max]
        weights = 0.5 * b_max * weights

        integral = 0.0
        for b, w in zip(b_scaled, weights):
            theta = nstp.scattering_angle(b, E_com, Z1, Z2, a_u, order=order)
            integral += w * np.sin(theta/2)**2 * b

        S_n = 2 * np.pi * gamma * E_lab * integral * 1e-16
        results[f"{order}"] = S_n

    return results

if __name__ == "__main__":
    Z1_H, Z2_Si = 1, 14
    M1_H, M2_Si = 1.00784, 28.0855
    Z1_Si, Z2_Au = 14, 79
    M1_Si, M2_Au = 28.0855, 196.96657
    order_list=[50,100, 200, 300, 400, 500]
    b_max_factor_list=[5,15,20,25,30,35,40,45,50]
    result = pd.DataFrame(columns=['collision', 'e_lab', 'b_max_factor', 'order', 's_n'])
    # Loop over different b_max factor
    for b_max_factor in b_max_factor_list:
        # Case 1: H --> Si, 10ev        
        results_H_Si = nuclear_stopping_ZBL_convergence_test(10, Z1_H, Z2_Si, M1_H, M2_Si,b_max_factor, order_list)
        for order, S_n in results_H_Si.items():
            result = result.append({'collision': 'H-Si', 'e_lab': 10, 'b_max_factor': b_max_factor, 'order': order, 's_n': S_n}, ignore_index=True)

        # Case 2: H --> Si,  5e6
        results_H_Si = nuclear_stopping_ZBL_convergence_test(5e6, Z1_H, Z2_Si, M1_H, M2_Si,b_max_factor, order_list)
        for order, S_n in results_H_Si.items():
            result = result.append({'collision': 'H-Si', 'e_lab': 5000000, 'b_max_factor': b_max_factor, 'order': order, 's_n': S_n}, ignore_index=True)

        # Case 3: Si --> Au; 10 ev
        results_Si_Au = nuclear_stopping_ZBL_convergence_test(10, Z1_Si, Z2_Au, M1_Si, M2_Au, b_max_factor, order_list)
        for order, S_n in results_Si_Au.items():
            result = result.append({'collision': 'Si-Au', 'e_lab': 10, 'b_max_factor': b_max_factor, 'order': order, 's_n': S_n}, ignore_index=True)
        
        # Case 4: Si --> Au; 5 Mev
        results_Si_Au = nuclear_stopping_ZBL_convergence_test(5e6, Z1_Si, Z2_Au, M1_Si, M2_Au, b_max_factor, order_list)
        for order, S_n in results_Si_Au.items():
            result = result.append({'collision': 'Si-Au', 'e_lab': 5000000, 'b_max_factor': b_max_factor, 'order': order, 's_n': S_n}, ignore_index=True)
    # Save the result to csv
    result.to_csv('../run/convergence_result.csv', index=False)  