#Plot relative error
# This script reads stopping power data for H-Si and Si-Au from CSV files, calculates the relative error between two methods (ZBL and universal), and plots the results.
# It uses the pandas library for data manipulation and matplotlib for plotting.
# The script also calculates the maximum relative error for both datasets and includes it in the plot legend.
# The plot is saved as a PNG file and displayed to the user.
# 

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    # H-Si
    df_H_Si = pd.read_csv('../run/H_Si_stopping_power.csv')
    h_si_energy = df_H_Si['Energy (eV)']
    h_si_sn_zbl = df_H_Si['S_n_H_Si_ZBL']
    h_si_sn_univ = df_H_Si['S_n_H_Si_univ']
    h_si_relative_error = np.abs(h_si_sn_zbl - h_si_sn_univ) / h_si_sn_zbl * 100
    # Si-Au
    df_Si_Au = pd.read_csv('../run/Si_Au_stopping_power.csv')
    si_au_energy = df_Si_Au['Energy (eV)']
    si_au_sn_zbl = df_Si_Au['S_n_Si_Au_ZBL']
    si_au_sn_univ = df_Si_Au['S_n_Si_Au_univ']
    si_au_relative_error = np.abs(si_au_sn_zbl - si_au_sn_univ) / si_au_sn_zbl * 100
    max_h_si_relative_error = np.max(h_si_relative_error)
    max_si_au_relative_error = np.max(si_au_relative_error)
    # Plotting
    plt.figure(figsize=(12, 8))
    plt.loglog(h_si_energy, h_si_relative_error, label=f'H-Si Relative Error,Max relative difference:{max_h_si_relative_error:0.3}%', color='blue')
    plt.loglog(si_au_energy, si_au_relative_error, label=f'Si-Au Relative Error,Max relative differenc:{max_si_au_relative_error:0.3}%', color='red')
    plt.xlabel('Energy (eV)', fontsize=12)
    plt.ylabel('Relative Error (%)', fontsize=12)
    plt.title('Relative Error in Nuclear Stopping Power', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.xlim(10, 5e6)
    plt.tight_layout()
    plt.savefig('../run/Relative_Error_Nuclear_Stopping_Power.png')
    plt.show()
    print("Relative error plot saved as Relative_Error_Nuclear_Stopping_Power.png")
