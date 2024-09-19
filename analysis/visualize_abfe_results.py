#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
===============================================================================
CopyRight (c) By freeenergylab.
@Description:
A tool can be used to visualize several replicates of ABFE results.

@Author: Pengfei Li
@Date: Aug. 18, 2024
===============================================================================
"""
import argparse
import os
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from copy import deepcopy
from scipy.constants import R, calorie

pd.options.display.max_rows=None
pd.options.display.float_format = "{:,.2f}".format

def f1(x): return x - 1.0

def f2(x): return x + 1.0

def f3(x): return x - 2.0

def f4(x): return x + 2.0

def beta(temperature):
    """Get beta (in mol/kJ) at this known temperature.
    """
    # set inverse temperature
    R_kcalmol = R/1000./calorie
    k_b = R_kcalmol
    beta = 1./(k_b*temperature) # kBT in kcal/mol

    return beta

class Analyzer(object):
    """Analyzer for visualizing these abfe results.

    Args:
        exp_file (str): the experimental result file.
    """
    def __init__(self, abfe_result_files, exp_file=None, method='MBAR',
                 shift=False, plot=False, figure_name='dG_comparison',
                 avg_type='direct', temperature=None):
        self.workdir = os.getcwd()
        if exp_file is not None:
            self.df_exp = pd.read_csv(exp_file, sep=',', index_col=0)
            self.df_exp.rename(columns={f'exp_dG': 'exp'}, inplace=True)
        self.df_abfe_list = []
        for i, file in enumerate(abfe_result_files):
            df = pd.read_csv(file, sep=',', index_col=0)
            df = pd.DataFrame(df, columns=[f'{method}_dG'])
            df.rename(columns={f'{method}_dG': f'run{i+1}'}, inplace=True)
            self.df_abfe_list.append(df)
        # merge several replicates of runs
        df_abfes = self.df_abfe_list[0]
        for n in range(1, len(self.df_abfe_list)):
            df_abfes = pd.merge(df_abfes, self.df_abfe_list[n], left_index=True, right_index=True)
        # merge the replicates of results and experimental results
        df_abfes = self.cal_avg_std(df=df_abfes, avg_type=avg_type, temperature=temperature)
        print('The mean of these <stdev> values is %.2f' % df_abfes.stdev.mean())
        print('The mean of these <spread> values is %.2f' % df_abfes.spread.mean())
        if exp_file is not None:
            self.df_abfes = pd.merge(df_abfes, self.df_exp, left_index=True, right_index=True)
            self.df_abfes['difference'] = df_abfes['average'] - self.df_abfes['exp']
            if shift:
                # do shift
                dG_shift = self.df_abfes.exp.mean() - self.df_abfes.average.mean()
                self.df_abfes.loc[:, "average"] = self.df_abfes["average"].apply(lambda x: x + dG_shift)
                self.df_abfes.loc[:, "difference"] = self.df_abfes["difference"].apply(lambda x: x + dG_shift)
                print('The mean of these <|difference|> values is %.2f' % self.df_abfes.difference.abs().mean())
                print('The following <average> values are shifted by %.2f kcal/mol.' % dG_shift)
            self.metrics = self.compute_metrics(df_input=self.df_abfes)
            self.metrics_std = self.bootstrap()
            # do plot
            if plot:
                self.do_plot(df_input=self.df_abfes, figure_name=figure_name)
        else:
            self.df_abfes = df_abfes

    def cal_avg_std(self, df, avg_type, temperature):
        """Calculate the average and stdev values for df."""
        df_copy = deepcopy(df)
        max_ = df_copy.max(axis=1, skipna=True)
        min_ = df_copy.min(axis=1, skipna=True)

        if avg_type == 'direct':
            df['average'] = df_copy.mean(axis=1, skipna=True)
            df['stdev'] = df_copy.std(axis=1, skipna=True)
        elif avg_type == 'boltzmann':
            df_new = df_copy.apply(lambda x: np.exp(-1.*beta(temperature)*x), axis=1)
            df['average'] = df_new.sum(axis=1, skipna=True)
            df.loc[:, "average"] = df["average"].apply(lambda x: -1./beta(temperature)*np.log(x))
            df['stdev'] = df_copy.std(axis=1, skipna=True)

        df['spread'] = max_ - min_

        return df

    def compute_metrics(self, df_input):
        """Compute the metrics between the predicted and experimental values.
        """
        df = pd.DataFrame(df_input, columns=['average', 'exp'])
        metrics = {}
        # pearson: pearson correlation coefficient
        metrics['r2_pearson'] = df['average'].corr(df['exp'], method='pearson')**2.0
        metrics['r_pearson'] = df['average'].corr(df['exp'], method='pearson')
        # spearman: spearman rank correlation coefficient
        metrics['rho_spearman'] = df['average'].corr(df['exp'], method='spearman')
        # kendall: kendall tau correlation coefficient
        metrics['tau_kendall'] = df['average'].corr(df['exp'], method='kendall')
        # root mean square error (deviation)
        metrics['rmse'] = ((df.average - df.exp)**2.0).mean()**0.5
        # mean signed error (deviation)
        metrics['mse'] = (df.average - df.exp).mean()
        # mean unsigned error (deviation)
        metrics['mue'] = (df.average - df.exp).abs().mean()

        return metrics

    def bootstrap(self, random_seed=2024, repeats=1000):
        """Bootstrap for calculating the stdev of metrics."""
        np.random.seed(random_seed)
        nmols = len(self.df_abfes.exp)
        mol_index = list(range(nmols))
        r2_pearson, r_pearson, rho_spearman, tau_kendall, rmse, mse, mue = [], [], [], [], [], [], []
        for repeat in range(repeats):
            index = np.random.choice(mol_index, nmols, replace=True)
            df_copy = deepcopy(self.df_abfes.iloc[index])
            metrics = self.compute_metrics(df_input=df_copy)
            r2_pearson.append(metrics['r2_pearson'])
            r_pearson.append(metrics['r_pearson'])
            rho_spearman.append(metrics['rho_spearman'])
            tau_kendall.append(metrics['tau_kendall'])
            rmse.append(metrics['rmse'])
            mse.append(metrics['mse'])
            mue.append(metrics['mue'])
        metrics_std = {}
        metrics_std['r2_pearson_std'] = np.std(r2_pearson)
        metrics_std['r_pearson_std'] = np.std(r_pearson)
        metrics_std['rho_spearman_std'] = np.std(rho_spearman)
        metrics_std['tau_kendall_std'] = np.std(tau_kendall)
        metrics_std['rmse_std'] = np.std(rmse)
        metrics_std['mse_std'] = np.std(mse)
        metrics_std['mue_std'] = np.std(mue)

        return metrics_std

    def do_plot(self, df_input, figure_name):
        """Plot the comparison between the predicted and experimental results."""
        df = df_input.dropna()

        exp_all_list = df['exp'].tolist()
        average_all_list = df['average'].tolist()
        std_all_list = df['stdev'].tolist()
    
        min_exp = min(exp_all_list)
        max_exp = max(exp_all_list)
        min_average = min(average_all_list)
        max_average = max(average_all_list)
        if (min_exp <= min_average):
            x_min = min_exp - 2.0
        else:
            x_min = min_average - 2.0
        if (max_exp >= max_average):
            x_max = max_exp + 2.0
        else:
            x_max = max_average + 2.0

        fig, ax = plt.subplots(figsize=(10,10))
        plt.errorbar(exp_all_list[:], average_all_list[:], yerr=std_all_list[:],
                     fmt='o', markerfacecolor='none', ecolor='black', ms=8,
                     color='blue', elinewidth=2, capsize=4)
        rows = [
            'MUE (kcal/mol)',
            'RMSE (kcal/mol)',
            'Pearson %s' % 'R',
            'Spearman %s' % (u'\u03C1'),
            'Kendall %s' % (u'\u03C4'),
            ]
        metric = [
            '%.2f' % self.metrics['mue'],
            '%.2f' % self.metrics['rmse'],
            '%.2f' % self.metrics['r_pearson'],
            '%.2f' % self.metrics['rho_spearman'],
            '%.2f' % self.metrics['tau_kendall'],
            ]
        metric_std = [
            '%.2f' % self.metrics_std['mue_std'],
            '%.2f' % self.metrics_std['rmse_std'],
            '%.2f' % self.metrics_std['r_pearson_std'],
            '%.2f' % self.metrics_std['rho_spearman_std'],
            '%.2f' % self.metrics_std['tau_kendall_std'],
            ]
        cell_text = []
        for i in range(len(rows)):
            cell_text.append([rows[i] + ': ' + str(metric[i]) + '+-' + str(metric_std[i])])
        plt.table(cellText=cell_text, edges='open', colLabels=None, bbox = [0.1, 0.7, 0.25, 0.25])
        # ax.spines['right'].set_color('none')
        # ax.spines['top'].set_color('none')
        ax.xaxis.set_minor_locator(MultipleLocator(1.0))
        ax.yaxis.set_minor_locator(MultipleLocator(1.0))
        ax.tick_params(axis='both', which='major', labelsize=18)
        t = np.arange(x_min, x_max, 1e-5)
        y1 = f1(t)
        y2 = f2(t)
        y3 = f3(t)
        y4 = f4(t)
        plt.fill_between(t, y3, y4, facecolor='lightgrey', interpolate=False)
        plt.fill_between(t, y1, y2, facecolor='grey', interpolate=False)
    
        plt.title('', loc='center')
        plt.xlabel(r'$\Delta G_{bind}^{experiment}$ [kcal/mol]', fontsize=20)
        plt.ylabel(r'$\Delta G_{bind}^{predicted}$ [kcal/mol]', fontsize=20)
        #plt.legend(loc='upper left', prop={'size':20}, shadow=False)
        plt.xlim(x_min, x_max)
        plt.ylim(x_min, x_max)
        plt.savefig(f'{figure_name}.pdf')

def _parse_args():
    parser = argparse.ArgumentParser(
        description='Parse arguments!',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument(
        '--abfe_result_files',
        default=['./run1/dG_results.csv', './run2/dG_results.csv', './run3/dG_results.csv'],
        type=str, nargs='*',
        help='The abfe result files from our FEP-SPell-ABFE workflow.'
        )
    parser.add_argument(
        '--exp_file', default=None,
        help='The experimental result file.'
        )
    parser.add_argument(
        '--method', type=str, default='MBAR',
        help='Which method results are visualized.'
        )
    parser.add_argument(
        '--shift', action='store_true', default=False,
        help='Whether to shift the abfe values or not.'
        )
    parser.add_argument(
        '--plot', action='store_true', default=False,
        help='Whether to plot the comparison between abfe and experimental results or not.'
        )
    parser.add_argument(
        '--figure_name', type=str, default='dG_comparison',
        help='The name for the saved figure.'
        )
    parser.add_argument(
        '--avg_type', type=str, default='direct', choices=['direct', 'boltzmann'],
        help='The type for calcuating the average and stdev.'
        )
    parser.add_argument(
        '--temperature', type=float, default=298.15,
        help='The temperature for calculating beta.'
        )

    return parser.parse_args()

if __name__ == '__main__':
    """Visualize several replicates of abfe results.
    """
    print(time.strftime("%c"))
    args = _parse_args()
    abfe_result_files = args.abfe_result_files
    exp_file = args.exp_file
    method = args.method
    shift = args.shift
    plot = args.plot
    figure_name = args.figure_name
    avg_type = args.avg_type
    temperature = args.temperature
    if plot:
        if exp_file is not None and os.path.exists(exp_file):
            print(f'Plot the comparison between abfe and experimental results.')
            print(exp_file)
        else:
            print(f'The exp_file does not exist.')
    else:
        print(f'Do not plot the comparison between abfe and experimental results.')
    print(f'These abfe result files are visualized:')
    for file in abfe_result_files:
        print(file)
        if not os.path.exists(file):
           print(f'The {file} does not exist.') 
    print(f'The {method} results are visualized.')
    # do analysis for visualzing abfe results
    analyzer = Analyzer(abfe_result_files=abfe_result_files, exp_file=exp_file,
                        method=method, shift=shift, plot=plot, figure_name=figure_name,
                        avg_type=avg_type, temperature=temperature)
    print(analyzer.df_abfes)
    if exp_file is not None and os.path.exists(exp_file):
        for name, metric in analyzer.metrics.items():
            print('%15s: %.2f +- %.2f' % (name, metric, analyzer.metrics_std[name+'_std']))