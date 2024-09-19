#!/bin/bash

module load anaconda3/FEP-SPell-ABFE 

# L99A plot
python visualize_abfe_results.py --abfe_result_files ../benchmark/L99A/run1/dG_results_sysmetry.csv ../benchmark/L99A/run2/dG_results_sysmetry.csv ../benchmark/L99A/run3/dG_results_sysmetry.csv --method "MBAR" --exp_file ../benchmark/L99A/exp/exp.csv --plot --figure_name "L99A_MBAR_Comparison" 1>L99A_MBAR_Comparison.dat
python visualize_abfe_results.py --abfe_result_files ../benchmark/L99A/run1/dG_results_sysmetry.csv ../benchmark/L99A/run2/dG_results_sysmetry.csv ../benchmark/L99A/run3/dG_results_sysmetry.csv --method "MBAR" --exp_file ../benchmark/L99A/exp/exp.csv --shift --plot --figure_name "L99A_MBAR_Comparison_shift" 1>L99A_MBAR_Comparison_shift.dat

# BRD4 plot
python visualize_abfe_results.py --abfe_result_files ../benchmark/BRD4/run1/dG_results.csv ../benchmark/BRD4/run2/dG_results.csv ../benchmark/BRD4/run3/dG_results.csv --method "MBAR" --exp_file ../benchmark/BRD4/exp/exp.csv --plot --figure_name "BRD4_MBAR_Comparison" 1>BRD4_MBAR_Comparison.dat
python visualize_abfe_results.py --abfe_result_files ../benchmark/BRD4/run1/dG_results.csv ../benchmark/BRD4/run2/dG_results.csv ../benchmark/BRD4/run3/dG_results.csv --method "MBAR" --exp_file ../benchmark/BRD4/exp/exp.csv --shift --plot --figure_name "BRD4_MBAR_Comparison_shift" 1>BRD4_MBAR_Comparison_shift.dat

# TYK2 plot
python visualize_abfe_results.py --abfe_result_files ../benchmark/TYK2/run1/dG_results.csv ../benchmark/TYK2/run2/dG_results.csv ../benchmark/TYK2/run3/dG_results.csv --method "MBAR" --exp_file ../benchmark/TYK2/exp/exp.csv --plot --figure_name "TYK2_MBAR_Comparison" 1>TYK2_MBAR_Comparison.dat
python visualize_abfe_results.py --abfe_result_files ../benchmark/TYK2/run1/dG_results.csv ../benchmark/TYK2/run2/dG_results.csv ../benchmark/TYK2/run3/dG_results.csv --method "MBAR" --exp_file ../benchmark/TYK2/exp/exp.csv --shift --plot --figure_name "TYK2_MBAR_Comparison_shift" 1>TYK2_MBAR_Comparison_shift.dat

# THROMBIN plot
python visualize_abfe_results.py --abfe_result_files ../benchmark/THROMBIN/run1/dG_results.csv ../benchmark/THROMBIN/run2/dG_results.csv ../benchmark/THROMBIN/run3/dG_results.csv --method "MBAR" --exp_file ../benchmark/THROMBIN/exp/exp.csv --plot --figure_name "THROMBIN_MBAR_Comparison" 1>THROMBIN_MBAR_Comparison.dat
python visualize_abfe_results.py --abfe_result_files ../benchmark/THROMBIN/run1/dG_results.csv ../benchmark/THROMBIN/run2/dG_results.csv ../benchmark/THROMBIN/run3/dG_results.csv --method "MBAR" --exp_file ../benchmark/THROMBIN/exp/exp.csv --shift --plot --figure_name "THROMBIN_MBAR_Comparison_shift" 1>THROMBIN_MBAR_Comparison_shift.dat
