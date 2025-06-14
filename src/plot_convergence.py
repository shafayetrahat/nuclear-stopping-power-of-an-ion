# This script is used to plot the convergence results of the nuclear stopping power for H-Si and Si-Au systems.
# It reads the data from a CSV file, processes it, and generates plots for different configurations.
# This script uses the pandas library for data manipulation and matplotlib for plotting.

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    # Load the data
    df = pd.read_csv('../run/convergence_result.csv')
    order_list=[50,100, 200, 300, 400, 500]
    b_max_factor_list=[5,15,20,25,30,35,40,45,50]
    # Sort the data by the collision
    df['collision'] = df['collision'].astype('category')
    df['collision'].cat.set_categories(['H-Si', 'Si-Au'], inplace=True)
    df.sort_values('collision', inplace=True)
    df_H_Si = df[df['collision'] == 'H-Si']
    df_Si_Au = df[df['collision'] == 'Si-Au']
    df_H_Si_elab_10 = df_H_Si[df_H_Si['e_lab'] == 10]
    df_H_Si_elab_5e6 = df_H_Si[df_H_Si['e_lab'] == 5000000]
    
    # Plot for H-Si
    plt.figure(figsize=(12, 8))
    for b_max_factor in b_max_factor_list:
        df_H_Si_b_max = df_H_Si_elab_10[df_H_Si_elab_10['b_max_factor'] == b_max_factor]
        df_H_Si_b_max.sort_values('order', inplace=True)
        plt.plot(df_H_Si_b_max['order'], df_H_Si_b_max['s_n'], label=f'H-Si, e_lab=10, b_max_factor={b_max_factor}')
    plt.xlabel('Order')
    plt.ylabel('Nuclear Stopping Power (eV/(atoms/cm^2))')
    plt.title('Nuclear Stopping Power for H-Si at e_lab=10 eV')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('../run/H_Si_e_lab_10.png')
    plt.clf()
    # Plot for e_lab=5e6
    for b_max_factor in b_max_factor_list:
        df_H_Si_b_max = df_H_Si_elab_5e6[df_H_Si_elab_5e6['b_max_factor'] == b_max_factor]
        df_H_Si_b_max.sort_values('order', inplace=True)
        plt.plot(df_H_Si_b_max['order'], df_H_Si_b_max['s_n'], label=f'H-Si, e_lab=5e6, b_max_factor={b_max_factor}')
    plt.xlabel('Order')
    plt.ylabel('Nuclear Stopping Power (eV/(atoms/cm^2))')
    plt.title('Nuclear Stopping Power for H-Si at e_lab=5e6 eV')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('../run/H_Si_e_lab_5e6.png')

######################################################
    df_Si_Au_elab_10 = df_Si_Au[df_Si_Au['e_lab'] == 10]
    df_Si_Au_elab_5e6 = df_Si_Au[df_Si_Au['e_lab'] == 5000000]
    plt.figure(figsize=(12, 8))
    for b_max_factor in b_max_factor_list:
        df_Si_Au_b_max = df_Si_Au_elab_10[df_Si_Au_elab_10['b_max_factor'] == b_max_factor]
        df_Si_Au_b_max.sort_values('order', inplace=True)
        plt.plot(df_Si_Au_b_max['order'], df_Si_Au_b_max['s_n'], label=f'Si-Au, e_lab=10, b_max_factor={b_max_factor}')
    plt.xlabel('Order')
    plt.ylabel('Nuclear Stopping Power (eV/(atoms/cm^2))')
    plt.title('Nuclear Stopping Power for Si-Au at e_lab=10 eV')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('../run/Si_Au_e_lab_10.png')
    plt.clf()
    # Plot for e_lab=5e6
    for b_max_factor in b_max_factor_list:
        df_Si_Au_b_max = df_Si_Au_elab_5e6[df_Si_Au_elab_5e6['b_max_factor'] == b_max_factor]
        df_Si_Au_b_max.sort_values('order', inplace=True)
        plt.plot(df_Si_Au_b_max['order'], df_Si_Au_b_max['s_n'], label=f'Si-Au, e_lab=5e6, b_max_factor={b_max_factor}')
    plt.xlabel('Order')
    plt.ylabel('Nuclear Stopping Power (eV/(atoms/cm^2))')
    plt.title('Nuclear Stopping Power for Si-Au at e_lab=5e6 eV')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('../run/Si_Au_e_lab_5e6.png')
    plt.clf()
    
################# Vary b_max_factor for fixed order ######################
    
    # Plot for H-Si
    plt.figure(figsize=(12, 8))
    for order in order_list:
        df_H_Si_order = df_H_Si_elab_10[df_H_Si_elab_10['order'] == order]
        df_H_Si_order.sort_values('b_max_factor', inplace=True)
        plt.plot(df_H_Si_order['b_max_factor'], df_H_Si_order['s_n'], label=f'H-Si, e_lab=10, order={order}')
    plt.xlabel('b_max_factor')
    plt.ylabel('Nuclear Stopping Power (eV/(atoms/cm^2))')
    plt.title('Nuclear Stopping Power for H-Si at e_lab=10 eV')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('../run/H_Si_e_lab_10_b_max_factor.png')
    plt.clf()
    # Plot for e_lab=5e6
    for order in order_list:
        df_H_Si_order = df_H_Si_elab_5e6[df_H_Si_elab_5e6['order'] == order]
        df_H_Si_order.sort_values('b_max_factor', inplace=True)
        plt.plot(df_H_Si_order['b_max_factor'], df_H_Si_order['s_n'], label=f'H-Si, e_lab=5e6, order={order}')
    plt.xlabel('b_max_factor')
    plt.ylabel('Nuclear Stopping Power (eV/(atoms/cm^2))')
    plt.title('Nuclear Stopping Power for H-Si at e_lab=5e6 eV')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('../run/H_Si_e_lab_5e6_b_max_factor.png')
    plt.clf()
    # Plot for Si-Au
    plt.figure(figsize=(12, 8))
    for order in order_list:
        df_Si_Au_order = df_Si_Au_elab_10[df_Si_Au_elab_10['order'] == order]
        df_Si_Au_order.sort_values('b_max_factor', inplace=True)
        plt.plot(df_Si_Au_order['b_max_factor'], df_Si_Au_order['s_n'], label=f'Si-Au, e_lab=10, order={order}')
    plt.xlabel('b_max_factor')
    plt.ylabel('Nuclear Stopping Power (eV/(atoms/cm^2))')
    plt.title('Nuclear Stopping Power for Si-Au at e_lab=10 eV')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('../run/Si_Au_e_lab_10_b_max_factor.png')
    plt.clf()
    # Plot for e_lab=5e6
    for order in order_list:
        df_Si_Au_order = df_Si_Au_elab_5e6[df_Si_Au_elab_5e6['order'] == order]
        df_Si_Au_order.sort_values('b_max_factor', inplace=True)
        plt.plot(df_Si_Au_order['b_max_factor'], df_Si_Au_order['s_n'], label=f'Si-Au, e_lab=5e6, order={order}')
    plt.xlabel('b_max_factor')
    plt.ylabel('Nuclear Stopping Power (eV/(atoms/cm^2))')
    plt.title('Nuclear Stopping Power for Si-Au at e_lab=5e6 eV')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('../run/Si_Au_e_lab_5e6_b_max_factor.png')
    plt.clf()