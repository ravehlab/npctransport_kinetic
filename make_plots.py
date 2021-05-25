import pickle
import sys
import numpy as np
import transport_simulation

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

def plot_MW_stats_list(stats_list_by_force, TSs_by_force):
    plot_from_sec=  0.1 # ts.bleach_start_time_sec + 1.0
    extras= [#'GTP_N',
            #'GDP_N',
            #'GTP_C',
            #'GDP_C',
            #'complexL_C',
            #'freeL_C',
            #'complexL_N',
            #'freeL_N'
    ]
    fig, ax_grid= plt.subplots(7+len(extras), 3,
                               figsize=(15,
                                        40.0 + 5.0*len(extras)),
                               sharex= False, sharey=False)
    n_NLS= len(stats_list_by_force[False])
    assert(n_NLS == len(stats_list_by_force[True]))
    ratios= np.ones(shape=(7+len(extras), n_NLS))
    ax_grid= ax_grid.transpose()
    for axes, is_force in zip(ax_grid[0:2,:], [False, True]):
        for i_NLS, stats in enumerate(stats_list_by_force[is_force]):
            ts= TSs_by_force[is_force][i_NLS]
            labels= ['L', 'U']
            x= stats['time_sec']
            ys={}
            ys[0]= stats['nuclear_importL_per_sec'] + stats['nuclear_importU_per_sec']
            ys[1]= stats['nuclear_exportL_per_sec'] + stats['nuclear_exportU_per_sec']
            ys[2]= get_N_C_ratio_stats(ts,
                                       stats,
                                       labels)
            ys[3]= get_compartment_concentration_stats(ts,
                                                      stats,
                                                      'C',
                                                      labels)
            ys[4]= get_compartment_concentration_stats(ts,
                                                      stats,
                                                      'N',
                                                      labels)
            ys[5]= ys[3] - ys[4]
            ys[6]= stats['complexL_NPC_N_import']+stats['complexL_NPC_C_import']+stats['complexL_NPC_N_export']+stats['complexL_NPC_C_export'] \
                + stats['complexU_NPC_C_import']+stats['complexU_NPC_C_import']+stats['complexU_NPC_N_export']+stats['complexU_NPC_C_export']
            for iextra, extra in enumerate(extras):
                ys[7+iextra]= stats[extra]
            plot_from_frame= int(plot_from_sec/ts.dt_sec)
            for i_row, ax in enumerate(axes):
                ax.plot(x[plot_from_frame:],
                        ys[i_row][plot_from_frame:],
                        label= free_to_complex_rates[i_NLS])
                ax.set_xlabel(r"time [$sec$]")
                if is_force:
                   ratios[i_row, i_NLS] *= ys[i_row][-1]
                else:
                   ratios[i_row, i_NLS] /= ys[i_row][-1]
            axes[0].set_ylabel(r"import rate [$sec^{-1}$]")
            axes[0].set_ylim([0.01,0.3])
            #axes[0].set_yscale('log')
            axes[1].set_ylabel(r"export rate [$sec^{-1}$]")
            axes[1].set_ylim([0.01,0.3])
            #axes[1].set_yscale('log')
            axes[2].set_ylabel("N/C ratio")
            axes[2].set_ylim([0,7.0])
            axes[3].set_ylabel(r"C [$M$]")
            axes[3].set_yscale('log')
            axes[4].set_ylabel(r"N [$M$]")
            axes[4].set_yscale('log')
            axes[5].set_ylabel(r"$\Delta$(C,N) [$M$]")
            axes[5].set_yscale('symlog', linthreshy=1e-9)
            axes[6].set_ylabel('NPC [nmol]')
            for iextra, extra in enumerate(extras):
                axes[7+iextra].set_ylabel(extra)
                axes[7+iextra].set_yscale('log')
            title= "30 kPa" if is_force else "5 kPa"
            axes[0].set_title(title)

    NLSs= [free_to_complex_rates[i_NLS] for i_NLS in range(ratios.shape[1])]
    ax_grid[2,0].set_title("Mechanosensitivity")
    for i_row, ax in enumerate(ax_grid[2,:]):
       ax.bar(range(len(NLSs)),
              ratios[i_row,:],
               width=0.8,
             tick_label= NLSs)
       ax.set_xlabel('NLS strength')

    handles, labels = ax_grid[0,0].get_legend_handles_labels()
    lh= fig.legend(handles, labels, loc='center left')
    lh.set_title('NLS strength')


def get_compartment_nmol_stats(ts, 
                               stats, 
                               compartment,
                              labels= ['L', 'U']):
    assert(compartment in ['N','C', 'NPC'])
    nframes= len(stats['time_sec'])
    nmol_stats= np.zeros(nframes)
    if compartment=='NPC':
        for label in labels:
            for side in ['N', 'C']:
                for source in ['import', 'export']:
                    tag= 'complex{}_NPC_{}_{}'.format(label,
                                                     side,
                                                     source)
                    nmol_stats = nmol_stats + stats[tag]
    else:
        for state in ['free','complex']:
            for label in labels:
                tag= '{}{}_{}'.format(state, 
                                      label, 
                                      compartment)
                nmol_stats = nmol_stats + stats[tag]
    return nmol_stats

def get_compartment_concentration_stats(ts, 
                                        stats, 
                                        compartment, 
                                       labels= ['L', 'U']):
    assert(compartment in ['N','C'])
    nmol_stats= get_compartment_nmol_stats(ts, 
                                           stats, 
                                           compartment,
                                           labels)
    is_nuclear= (compartment=='N')
    volume_L= (ts.get_v_N_L() if is_nuclear else ts.get_v_C_L())
    return (nmol_stats/transport_simulation.N_A)/volume_L

def get_N_C_ratio_stats(ts, 
                        stats,
                        labels= ['L','U']):
    EPSILON= 1E-12
    c_N_stats= get_compartment_concentration_stats(ts, 
                                                   stats, 
                                                   'N')
    c_C_stats= get_compartment_concentration_stats(ts, 
                                                   stats, 
                                                   'C')
    return c_N_stats/c_C_stats


def make_plot(stats_by_force, figname, free_to_complex_rates):
    plot_MW_stats_list(*stats_by_force, free_to_complex_rates)
    plt.savefig(figname)
    
