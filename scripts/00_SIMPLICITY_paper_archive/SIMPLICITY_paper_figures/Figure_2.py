# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn           as sns
import matplotlib.gridspec as gridspec

import simplicity.dir_manager      as dm
import simplicity.output_manager   as om

import simplicity.tuning.evolutionary_rate as er
import simplicity.plots_manager as pm
from simplicity.intra_host_model import Host

pm.apply_plos_rcparams()

def plot_intra_host(ax):
    """
    Plot intra-host model solutions over time for all initial states.

    Displays curves for the probability of being infectious at time t for
    all possible initial states (0 to 20), with a continuous color mapping.

    """
    ax.text(-0.1, 1.05, "A", transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top', ha='left')
    ih_model = Host(tau_1=2.86, tau_2=3.91, tau_3=7.5, tau_4=8)
    time = 50
    step = 0.1
    
    t = np.arange(0, time, step)
    # Plot probability curves per state
    
    p_inf, p_det, p_rec = ih_model.data_plot_ih_solution(0,time,step)
    
    ax.plot(t, p_inf, color='red', label='Probability of being infectious')
    ax.plot(t, p_det, color='blue', label='Probability of being detectable')
    ax.plot(t, p_rec, color='black', label='Probability of being recovered')

    # Axis labels
    ax.set_xlabel("Time (d) ")
    ax.set_ylabel("Probability of being infected after t days")
    ax.set_xlim(0, time)
    ax.set_ylim(0, 1.3)
    ax.legend(loc = 'upper left')

    pm.apply_standard_axis_style(ax)


def plot_figure_tempest_regression(experiment_name, ax, out_sim_dir_i=9):
    
    ax.text(-0.1, 1.05, "B", transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top', ha='left')
    # parameter = 'nucleotide_substitution_rate'
    simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
    simulation_output_dir = simulation_output_dirs[out_sim_dir_i]
    # param = sm.get_parameter_value_from_simulation_output_dir(simulation_output_dir, parameter)
    
    sequencing_data_df = om.create_combined_sequencing_df(simulation_output_dir, 
                                                          min_seq_number=30,
                                                          min_sim_lenght=100)
    

    fitted_tempest_regression = er.tempest_regression(sequencing_data_df)
    
    pm.plot_tempest_regression(sequencing_data_df,
                               fitted_tempest_regression,
                               ax)
    pm.apply_standard_axis_style(ax)

def plot_OSR_fit_figure(experiment_name, ax):
    ''' plot fit of nucleotide substitution rate / observed substitution rates curve
    '''
    ax.text(-0.1, 1.05, "C", transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top', ha='left')
    parameter = 'nucleotide_substitution_rate'
    model_type = 'exp'
    min_seq_number=30
    min_sim_lenght=100
    weights = None
    
    line_color = 'black' #'#DE8F05' # orange
    scatter_color = '#DE8F05'# orange
    combined_OSR_marker_color = '#E64B9D' # pink
    scatter_color_2 = '#0173B2' # blue '#029E73' # green
    
    # build the dataframe needed for the fit
    om.write_combined_OSR_vs_parameter_csv(experiment_name, 
                                            parameter, 
                                            min_seq_number,
                                            min_sim_lenght)
    
    om.write_OSR_vs_parameter_csv(experiment_name, 
                                    parameter, 
                                    min_seq_number,
                                    min_sim_lenght)
    
    # import the df needed for the fit
    OSR_single_sim_data = om.read_OSR_vs_parameter_csv(experiment_name, 
                                         parameter,
                                         min_seq_number,
                                         min_sim_lenght)
    
    fit_result = er.fit_observed_substitution_rate_regressor(experiment_name,
                                                             OSR_single_sim_data, 
                                                             model_type, weights)

    # import combined regression data
    OSR_combined_sim_data = om.read_combined_OSR_vs_parameter_csv(experiment_name,
                                                         parameter,
                                                         min_seq_number,
                                                         min_sim_lenght)
   
    # Group by nucleotide_substitution_rate and compute mean and standard deviation for OSR
    OSR_mean_data = om.get_mean_std_OSR(experiment_name,
                                       parameter,
                                       min_seq_number,
                                       min_sim_lenght)
    
    
    x_data = OSR_single_sim_data[parameter]
    
    # shaded area for real SARS-CoV-2 rates range
    ax.axhspan(1e-4, 1e-3, color='lightgray', alpha=0.2, zorder=-1, label= 'SARS-CoV-2 OSR lit. data range')
   
    # scatterplot Estimated OSR - single simulation
    sns.scatterplot(x=parameter, y='observed_substitution_rate', 
                    label='Estimated OSR - single simulation', data=OSR_single_sim_data,
                    color=scatter_color, alpha=0.5, ax=ax,
                    zorder=0)
    # plot fitted curve
    sns.lineplot(x=x_data, y=fit_result.best_fit, 
                 label=f'Fitted {model_type} curve', 
                 color=line_color, linewidth=2, ax=ax,
                 zorder=1)
    
    # scatterplot combined regression points (as comparison)
    sns.scatterplot(x=parameter, y='observed_substitution_rate', marker='X',
        label='Combined tempest regression estimate of OSR', data=OSR_combined_sim_data,
        color=combined_OSR_marker_color,alpha=1, ax=ax,
        zorder=2)
    
    # plot mean of observed_substitution_rate
    sns.scatterplot(x=parameter, y='mean', marker = 'X',
                    label='Mean of estimated OSR (single simulations)', data=OSR_mean_data,
                    color=scatter_color_2, alpha=1, ax=ax,
                    zorder=3)
    
    # SHOW VALUES USED IN SIMULATIONS
    A = fit_result.params['A'].value
    B = fit_result.params['B'].value
    C = fit_result.params['C'].value

    used_nsr = 0.00008759

    # compute the corresponding OSR 
    OSR_target = A * (used_nsr ** B) + C

    # 4) draw a hollow red circle at (x=target_nsr, y=OSR_target):
    ax.scatter(
        used_nsr,
        OSR_target,
        s=150,                    
        facecolors='none',
        edgecolors='red',
        linewidth=2,
        label= "NSR/OSR used for experiments",
        zorder=4
    )
    
    # Set axis (log log scale) -----------------------------------------------
    ax.set_xlabel('Nucleotide substitution rate (ssy)')
    ax.set_ylabel('Observed substitution rate (ssy)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(x_data.min(),x_data.max()*1.1)
    ax.legend()
    
    pm.apply_standard_axis_style(ax)
    

def plot_figure_2(experiment_name):
    fig = plt.figure(figsize=(20, 15))
    gs = gridspec.GridSpec(2, 2, figure=fig)

    # Subplot 1: Intra-host model 
    ax1 = fig.add_subplot(gs[0, 0])
    plot_intra_host(ax1)

    # Subplot 2: Tempest regression 
    ax2 = fig.add_subplot(gs[0, 1])
    plot_figure_tempest_regression(experiment_name, ax2)

    # Subplot 3: OSR fit 
    ax3 = fig.add_subplot(gs[1, 0:2])
    plot_OSR_fit_figure(experiment_name, ax3)
     
    fig.subplots_adjust(top=0.95, bottom=0.05)
    
    output_path = os.path.join("Data", f"Figure_2_{experiment_name}.tiff")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    output_path = os.path.join("Data", "Figure_2.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to: {output_path}")
    plt.close(fig)


def main(): 
    
    plot_figure_2(experiment_name="OSR_fit")
    
if __name__ == '__main__':
    main() 
    
    