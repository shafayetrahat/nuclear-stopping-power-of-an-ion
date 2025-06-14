"""
@Author: Md Shafayet Islam
@Date: 2025-05-18
@Description:
This program calculates the nuclear stopping power of ions in matter using two different models:
    1. The Ziegler-Biersack-Littmark (ZBL) model
    2. The universal model.
It computes the stopping power for two cases:
    1. Hydrogen ions (H) in Silicon (Si)
    2. Silicon ions (Si) in Gold (Au).
The program uses numerical integration to evaluate the scattering angle and the stopping power.It 
also plots the results and calculates the relative differences between the two models.
@Usage:
    python main.py
@Output:
    - A plot comparing the nuclear stopping power for both cases.
    - The maximum relative differences between the two models for each case.
@Dependencies:
    - numpy
    - scipy
    - matplotlib
    - Python 3.x
"""

# Import necessary libraries
import numpy as np
from scipy.optimize import bisect
import matplotlib.pyplot as plt
import pandas as pd
# Constants
e_sq = 14.3996  # e^2/(4 pi epsilon_0) in eV·Angstrom
alpha = np.array([0.1818, 0.5099, 0.2802, 0.02817])
beta = np.array([3.2, 0.9423, 0.4028, 0.2016])

def screening_function(x):
    """
    Calculate the screening function for the ZBL model.
    Args:
        x (float): The ratio of distance to screening length.
    Returns:
        float: The screening function value.
    """
    return np.sum(alpha * np.exp(-beta * x))

def potential(r, Z1, Z2, a_u):
    """
    Calculate the potential energy between two ions.
    Args:
        r (float): Distance between ions in Angstrom.
        Z1 (int): Atomic number of ion 1.
        Z2 (int): Atomic number of ion 2.
        a_u (float): Screening length in Angstrom.
    Returns:
        float: The potential energy in eV.
    """
    if r == 0:
        return np.inf
    return (Z1 * Z2 * e_sq / r) * screening_function(r / a_u)

def g_tilde(r, b, E_com, Z1, Z2, a_u):
    """
    Calculate the function g_tilde for the given parameters.
    Args:
        r (float): Distance between ions in Angstrom.
        b (float): Impact parameter in Angstrom.
        E_com (float): Center of mass energy in eV.
        Z1 (int): Atomic number of ion 1.
        Z2 (int): Atomic number of ion 2.
        a_u (float): Screening length in Angstrom.
    Returns:
        float: The value of g_tilde.
    """
    return 1 - (b/r)**2 - potential(r, Z1, Z2, a_u)/E_com

def find_r_min(b, E_com, Z1, Z2, a_u):
    """
    Find the minimum distance of approach for the given parameters.
    Args:
        b (float): Impact parameter in Angstrom.
        E_com (float): Center of mass energy in eV.
        Z1 (int): Atomic number of ion 1.
        Z2 (int): Atomic number of ion 2.
        a_u (float): Screening length in Angstrom.
    Returns:
        float: The minimum distance of approach in Angstrom.
    """
    try:
        r_min = bisect(g_tilde, b, 100*b, args=(b, E_com, Z1, Z2, a_u), rtol=1e-6)
        return r_min
    except ValueError:
        return np.nan  # Fallback for grazing collisions

def integrand(u, b, r_min, E_com, Z1, Z2, a_u):
    """
    Integrand for the scattering angle calculation.
    Args:
        u (float): Variable of integration.
        b (float): Impact parameter in Angstrom.
        r_min (float): Minimum distance of approach in Angstrom.
        E_com (float): Center of mass energy in eV.
        Z1 (int): Atomic number of ion 1.
        Z2 (int): Atomic number of ion 2.
        a_u (float): Screening length in Angstrom.
    Returns:
        float: The value of the integrand.
    """
    u_trans = u  # Using quad directly on [0,1]
    if u_trans < 1e-6:
        return 0
    r = r_min / (1 - u_trans**2)
    term = b**2 * (2 - u_trans**2) + (r_min**2 / (u_trans**2 * E_com)) * (
        potential(r_min, Z1, Z2, a_u) - potential(r, Z1, Z2, a_u))
    return 1/np.sqrt(term) if term > 0 else 0

def scattering_angle(b, E_com, Z1, Z2, a_u, order):
    """
    Calculate the scattering angle using Gauss-Legendre quadrature.
    Args:
        b (float): Impact parameter in Angstrom.
        E_com (float): Center of mass energy in eV.
        Z1 (int): Atomic number of ion 1.
        Z2 (int): Atomic number of ion 2.
        a_u (float): Screening length in Angstrom.
        order (int): Order of the Gauss-Legendre quadrature.
    Returns:
        float: The scattering angle in radians.
    """
    r_min = find_r_min(b, E_com, Z1, Z2, a_u)
    if np.isnan(r_min):
        return 0.0
    
    nodes, weights = np.polynomial.legendre.leggauss(order)
    u_nodes = 0.5 * (nodes + 1)
    u_weights = 0.5 * weights
    
    # Vectorized computation
    mask = u_nodes >= 1e-6  # Filter out near-zero values
    u_nodes = u_nodes[mask]
    u_weights = u_weights[mask]
    
    # Vectorized calculations
    r_vals = r_min / (1 - u_nodes**2)
    potential_diff = potential(r_min, Z1, Z2, a_u) - np.vectorize(potential)(r_vals, Z1, Z2, a_u)
    terms = b**2 * (2 - u_nodes**2) + (r_min**2 / (u_nodes**2 * E_com)) * potential_diff
    
    # Only consider positive terms
    valid_terms = terms > 0
    integral = np.sum(u_weights[valid_terms] / np.sqrt(terms[valid_terms]))
    
    return np.pi - 4 * b * integral

def nuclear_stopping_ZBL(E_lab, Z1, Z2, M1, M2, order=200, b_max_factor=25):
    """
    Calculate the nuclear stopping power using the ZBL model.
    Args:
        E_lab (float): Energy in eV.
        Z1 (int): Atomic number of ion 1.
        Z2 (int): Atomic number of ion 2.
        M1 (float): Mass of ion 1 in atomic mass units.
        M2 (float): Mass of ion 2 in atomic mass units.
        order (int): Order of the Gauss-Legendre quadrature.
        b_max_factor (float): Multiplier for screening length (a_u) to set b_max.
    Returns:
        float: The nuclear stopping power in eV/(atoms/cm^2).
    """
    E_com = (M2 / (M1 + M2)) * E_lab
    a_u = 0.46848 / (Z1**0.23 + Z2**0.23)
    gamma = 4 * M1 * M2 / (M1 + M2)**2
    
    # Get quadrature points
    b_nodes, b_weights = np.polynomial.legendre.leggauss(order)
    b_max = b_max_factor * a_u
    b_scaled = 0.5 * b_max * (b_nodes + 1)
    b_weights = 0.5 * b_max * b_weights
    
    # Vectorized scattering angle calculation
    thetas = np.array([scattering_angle(b, E_com, Z1, Z2, a_u, order) for b in b_scaled])
    
    # Vectorized integral computation
    integrand = b_weights * np.sin(thetas/2)**2 * b_scaled
    integral = np.sum(integrand)
    
    return 2 * np.pi * gamma * E_lab * integral * 1e-16

def nuclear_stopping_universal(E_lab, Z1, Z2, M1, M2):
    """
    Calculate the nuclear stopping power using the universal model.
    Args:
        E_lab (float): Energy in eV.
        Z1 (int): Atomic number of ion 1.
        Z2 (int): Atomic number of ion 2.
        M1 (float): Mass of ion 1 in atomic mass units.
        M2 (float): Mass of ion 2 in atomic mass units.
    Returns:
        float: The nuclear stopping power in eV/(atoms/cm^2).
    """
    E_com = (M2 / (M1 + M2)) * E_lab/1e3  # keV
    denominator = Z1 * Z2 * (Z1**0.23 + Z2**0.23)
    epsilon = 32.53 * E_com / denominator
    
    if epsilon <= 30:
        numerator = np.log(1 + 1.138*epsilon)
        denominator_sn = 2*(epsilon + 0.01321*epsilon**0.21226 + 0.19593*epsilon**0.5)
    else:
        numerator = np.log(epsilon)
        denominator_sn = 2*epsilon
    
    s_n = numerator / denominator_sn
    prefactor = 8.462e-15 * Z1 * Z2 * M1 / ((M1 + M2) * (Z1**0.23 + Z2**0.23))
    return prefactor * s_n

if __name__=="__main__":
    energy_points = 30
    # Energy range (eV)
    energies = np.logspace(1, np.log10(5e6), energy_points)

    # Case 1: H --> Si
    order,b_max_factor = 400, 15  # Converged value from the test
    Z1_H, Z2_Si = 1, 14
    M1_H, M2_Si = 1.00784, 28.0855
    S_n_H_Si_ZBL = np.array([nuclear_stopping_ZBL(E, Z1_H, Z2_Si, M1_H, M2_Si,order,b_max_factor) for E in energies])
    S_n_H_Si_univ = np.array([nuclear_stopping_universal(E, Z1_H, Z2_Si, M1_H, M2_Si) for E in energies])

    # Case 2: Si --> Au
    order,b_max_factor = 200, 30  # Converged value from the test
    Z1_Si, Z2_Au = 14, 79
    M1_Si, M2_Au = 28.0855, 196.96657
    S_n_Si_Au_ZBL = np.array([nuclear_stopping_ZBL(E, Z1_Si, Z2_Au, M1_Si, M2_Au,order,b_max_factor) for E in energies])
    S_n_Si_Au_univ = np.array([nuclear_stopping_universal(E, Z1_Si, Z2_Au, M1_Si, M2_Au) for E in energies])
    
    # Convert to dataframe
    data_H_Si = {
        'Energy (eV)': energies,
        'S_n_H_Si_ZBL': S_n_H_Si_ZBL,
        'S_n_H_Si_univ': S_n_H_Si_univ,
    }
    data_Si_Au = {
        'Energy (eV)': energies,
        'S_n_Si_Au_ZBL': S_n_Si_Au_ZBL,
        'S_n_Si_Au_univ': S_n_Si_Au_univ
    }
    df_H_Si = pd.DataFrame(data_H_Si)
    df_Si_Au = pd.DataFrame(data_Si_Au)
    df_H_Si.to_csv('../run/H_Si_stopping_power.csv', index=False)
    df_Si_Au.to_csv('../run/Si_Au_stopping_power.csv', index=False)
    print("Data saved to H_Si_stopping_power.csv and Si_Au_stopping_power.csv")
    # Ploting
    plt.figure(figsize=(12, 8))
    plt.loglog(energies, S_n_H_Si_ZBL, label='ZBL: H --> Si', color='blue')
    plt.loglog(energies, S_n_H_Si_univ, '--', label='Universal: H --> Si', color='blue')
    plt.loglog(energies, S_n_Si_Au_ZBL, label='ZBL: Si --> Au', color='red')
    plt.loglog(energies, S_n_Si_Au_univ, '--', label='Universal: Si --> Au', color='red')

    plt.xlabel('Energy (eV)', fontsize=12)
    plt.ylabel('Nuclear Stopping Power (eV/(atoms/cm^2))', fontsize=12)
    plt.title('Nuclear Stopping Power Comparison', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.xlim(10, 5e6)
    plt.tight_layout()
    plt.savefig('../run/Nuclear_Stopping_Power_Comparison.png')

    # Relative differences
    rel_diff_H_Si = np.abs(S_n_H_Si_ZBL - S_n_H_Si_univ) / S_n_H_Si_ZBL * 100
    rel_diff_Si_Au = np.abs(S_n_Si_Au_ZBL - S_n_Si_Au_univ) / S_n_Si_Au_ZBL * 100

    print(f"Max relative difference for H-->Si: {np.max(rel_diff_H_Si):.2f}%")
    print(f"Max relative difference for Si-->Au: {np.max(rel_diff_Si_Au):.2f}%")
