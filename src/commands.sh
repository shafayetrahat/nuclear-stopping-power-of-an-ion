#!/bin/bash
# This script runs the Python scripts to generate the plots for the paper.
# Dependencies:
# - Python 3.x
# - matplotlib
# - numpy
# - pandas
# - scipy
#!!!!!!!!!!!!!! CAUTION. Need high time to compute. !!!!!!!!!!!!!!
python3 convergence_test.py
python plot_convergence.py

##!!!!!!!!!!!!!! CAUTION. Need high time to compute. You can decrease the 'energy_points' in nuclear_stopping.py 
#                for more faster calculation                                                                    !!!!!!!!!!!!!!
python3 nuclear_stopping.py
python3 plot_relative_error.py