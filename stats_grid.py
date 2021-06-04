import multiprocessing
import numpy as np
import argparse
import map_param_grid
import matplotlib as mpl
import matplotlib.pyplot as plt
import transport_simulation 
from transport_simulation import TransportSimulation

main_description = "Create stats grids, with parameters fraction_complex_NPC_traverse_per_sec, rate_free_to_complex_per_sec, max_passive_diffusion_rate_nmol_per_sec_per_M, over a chosen range."
epilogue = "For example, running 'python stats_grid.py -c 50e-9 -nx 20 -ny 20 -n 0 3 -p 0.001 0.01 -np 12"\
           " -o example_fig' will create 12 figures (-np) with resolution 20x20 (-nx, -ny) from 20x20=400"\
           " simulations, with fraction_complex_NPC_traverse_per_sec ranging between log(0) and log(3)"\
           " (-n 0 3), and cargo concentration of 50e-9 M, each figure will have increasing passive diffusion"\
           " rates, in the range 0.001-0.01 (-p), and the figures will be saved as 'example_fig_p_{passive}.png"\
           " (-o example_fig )for each passive value in that range."


def get_param_range_traverse_kon(nx, ny, npc_traverse_range=(0,3), k_on_range=(-2,1)):
    param_range= {}
    print(f"nx={nx} ny={ny}")
    param_range['tag_x']= "fraction_complex_NPC_traverse_per_sec"
    param_range['range_x']= np.logspace(*npc_traverse_range, nx) 
    param_range['pretty_x']= r"rate NPC traverse [$sec^{-1}$]"
    param_range['tag_y']= "rate_free_to_complex_per_sec"
    param_range['range_y']= np.logspace(*k_on_range, ny)
    param_range['pretty_y']= r"NTR $k_{on}$ [$sec^{-1}$]"
    return param_range
def get_transport_simulation_by_passive(passive_nuclear_molar_rate_per_sec,
                                        Ran_cell_M = 20.0e-6,
                                        **kwargs):    
    v_N_L=627e-15
    v_C_L=2194e-15    
    ts= transport_simulation.TransportSimulation(v_N_L= v_N_L,
                                                v_C_L= v_C_L)
    ts.set_time_step(0.1e-3)
    ts.set_NPC_dock_sites(n_NPCs= 2000, 
                        n_dock_sites_per_NPC= 500)
    ts.set_passive_nuclear_molar_rate_per_sec(passive_nuclear_molar_rate_per_sec) #get_passive_export_rate_per_sec(27,1))
    ts.fraction_complex_NPC_to_free_N_per_M_GTP_per_sec = 1.0e+6 # TODO: this is doubled relative to complex_N to free_N
    ts.fraction_complex_N_to_free_N_per_M_GTP_per_sec = 1.0e+6
    ts.rate_complex_to_NPC_per_free_site_per_sec_per_M= 50e+6
    ts.fraction_complex_NPC_to_complex_N_C_per_sec= 3000.0 # Leakage parameter
    ts.rate_GDP_N_to_GTP_N_per_sec= 1000.0
    ts.rate_GTP_N_to_GDP_N_per_sec= 0.2
    ts.rate_GTP_C_to_GDP_C_per_sec= 500.0
    ts.rate_GTP_N_to_GTP_C_per_sec = 0.5
    ts.rate_GDP_C_to_GDP_N_per_sec = 1.0
    ts.rate_GDP_N_to_GDP_C_per_sec = 1.0 
    ts.rate_complex_to_free_per_sec = 0.05
    ts.rate_free_to_complex_per_sec = 0.01 # SCAN
    ts.fraction_complex_NPC_traverse_per_sec=4000 # SCAN
    ts.set_params(**kwargs) # override defaults
    ts.set_RAN_distribution(Ran_cell_M= Ran_cell_M, # total physiological concentration of Ran # TODO: check in the literature 
                                  parts_GTP_N=1000,
                                  parts_GTP_C=1,
                                  parts_GDP_N=1,
                                  parts_GDP_C=1000)
    return ts

def plot_stats_grids(stats_grids, transport_simulation, param_range,
                        NC_min=1.0,
                        NC_max= 20.0,
                        vmax_import_export=10.0):
    fig, axes= plt.subplots(3,3, figsize=(14, 10), sharex=True, sharey=True)
    # N/C
    plt.sca(axes[0,0])
    map_param_grid.plot_NC_ratios(param_range, 
                                  stats_grids, 
                                  transport_simulation, 
                                  vmin= NC_min,
                                  vmax= NC_max,
#                                  levels= np.linspace(NC_min, NC_max, 21)
                                  levels= np.logspace(np.log2(NC_min),np.log2(NC_max),21, base=2.0),
                                  locator= mpl.ticker.LogLocator(base=2.0)
                                 )
    # Bound fraction
    map_param_grid.plot_bound_fraction(param_range, stats_grids, 
                                       'N', ax=axes[0,1])
    map_param_grid.plot_bound_fraction(param_range, stats_grids, 
                                       'C', ax=axes[0,2])
    # Import/export
    import_export_locator= mpl.ticker.LogLocator(subs=[1.0, 5.0])
    map_param_grid.plot_import_export(param_range,
                                      stats_grids,
                                      axes= [axes[1,0], axes[1,1]],
                                   vmin= 0.01,
                                   vmax= vmax_import_export,
#                    levels=np.linspace(0.0,vmax_import_export,20),
                    levels= np.logspace(np.log10(1e-2),np.log10(vmax_import_export),21),
                                      locator= import_export_locator, 
                                     extend='both')
    plt.sca(axes[1,2])
    ratios_import_export= map_param_grid.get_import_export_ratios(stats_grids)
    map_param_grid.plot_param_grid(param_range, 
                    ratios_import_export,
                    Z_label= 'import:export',
                    vmin= NC_min,
                    vmax= NC_max,
#                   levels= np.linspace(NC_min, NC_max, 21)
                    levels= np.logspace(np.log2(NC_min),np.log2(NC_max),21, base=2.0),
                    locator= mpl.ticker.LogLocator(base=2.0) ,                                  
                    extend='both'
                                  )
    
    plt.sca(axes[2,0])
    GTP_ratio= stats_grids['GTP_N']/stats_grids['GTP_C']
    map_param_grid.plot_param_grid(param_range, 
                    GTP_ratio,
                    Z_label= "GTP N:C",
                    vmin=0.0, 
                     vmax=2000.0, 
                     levels=np.linspace(0,2000,21), 
                     extend='both')
    plt.sca(axes[2,1])
    GDP_ratio= stats_grids['GDP_C']/stats_grids['GDP_N']
    map_param_grid.plot_param_grid(param_range, 
                    GDP_ratio,
                    Z_label= f"GDP C:N",
                    vmin=0.0, 
#                     vmax=vmax, 
#                     levels=np.linspace(0,vmax,11), 
                     extend='both')
    plt.sca(axes[2,2])
    Ran_ratio= (stats_grids['GTP_N']+stats_grids['GDP_N'])/(stats_grids['GTP_C']+stats_grids['GDP_C'])
    map_param_grid.plot_param_grid(param_range, 
                    Ran_ratio,
                    Z_label= f"Ran N:C",
                    vmin=0.0, 
                     vmax=20, 
                    levels=np.linspace(0,20,21), 
                     extend='both')

def main_test(output, passive_range, npc_traverse_range, k_on_range, nx=20, ny=20, n_passive=10, pkl=None):
    print(output, passive_range, npc_traverse_range, k_on_range, nx, ny, n_passive, pkl)

def transport_simulation_generator(passive, Ran_cell_M, c_M, **kwargs):
    print(f"Ran: {Ran_cell_M:.6f} M")
    return get_transport_simulation_by_passive(passive_nuclear_molar_rate_per_sec= passive, 
                                               Ran_cell_M = Ran_cell_M, init_cargo_cytoplasm_M=c_M,
                                               **kwargs)

    
def main(output,
         passive_range,
         npc_traverse_range,
         k_on_range,
         nx=20, ny=20,
         n_passive=10,
         c_M=50e-6,
         pkl=None):
    test_ts= get_transport_simulation_by_passive(0.02, False)
    print(test_ts.max_passive_diffusion_rate_nmol_per_sec_per_M)
    Ran_cell_M = 80.0e-6
    param_range= get_param_range_traverse_kon(nx, ny, npc_traverse_range, k_on_range)
    print(param_range)
    n_processors= multiprocessing.cpu_count()
    stats_grids_traverse_by_passive_force= {} # 2D maps of statistics for different passive diffusion params
    ts_traverse_by_passive_force= {} #  transport simulaiton object used for each
    print("*** Starting multiprocess run ***")    
    
    for passive in np.logspace(*np.log10(passive_range), n_passive): #0.01,0.09,6):
            key= passive
            if key in stats_grids_traverse_by_passive_force:
                continue
            tsg_params = {"passive":passive,
                          "Ran_cell_M":Ran_cell_M,
                          "c_M": c_M}
            stats_grids_traverse_by_passive_force[key], \
            ts_traverse_by_passive_force[key] = \
                map_param_grid.map_param_grid_parallel\
                ( param_range,
                  equilibration_time_sec= 100.0,
                  n_processors= n_processors,
                  transport_simulation_generator= transport_simulation_generator,
                  transport_simulation_generator_params= tsg_params)
            print(f"passive rate {passive}")
            plot_stats_grids(stats_grids_traverse_by_passive_force[key],
                        ts_traverse_by_passive_force[key],
                        param_range,
                        vmax_import_export= 10.0,
                        NC_max=30.0,
                        NC_min=1.0)
            filename = f"{output}_p_{key}.png"
            plt.savefig(filename)
            print(f"Saved to {filename}")
            
            
    print("*** Finished multiprocess run ***")

    if pkl is not None:
    # Pickle results
        with open(pkl, "wb") as F:
            pickle.dump([stats_grids_traverse_by_passive_force, ts_traverse_by_passive_force], F)
            
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=main_description, epilog=epilogue)
    parser.add_argument("-np", "--n-passive", type=int, default=10, \
                        help="the number of values in the range --pasive-range to run")
    parser.add_argument("-p", "--passive-range", metavar=("MIN_PASSIVE","MAX_PASSIVE"), \
                        help="The range of values for passive_nuclear_molar_rate_per_sec", \
                        type=float, nargs=2, default=(0.01, 0.1))
    parser.add_argument("-n", "--npc-traverse-range", metavar=("MIN_RATE","MAX_RATE"), \
                        help="The range of values for fraction_complex_NPC_traverse_per_sec in log scale", \
                        type=float, nargs=2, default=(0., 3.0))
    parser.add_argument("--k-on-range", metavar=("MIN_K_ON", "MAX_K_ON"), \
                        help="The range of values for free_to_complex_per_sec in log scale", \
                        type=float, nargs=2, default=(-2., 1.0))
    parser.add_argument("-nx", type=int, default=20, \
                        help="Number of NPC_traverse_per_sec values to run in the range defined "\
                             "by --npc_traverse_range")
    parser.add_argument("-ny", type=int, default=20,
                        help="Number of rate_free_to_complex_per_sec values to run in the range "\
                             "defined by --k-on-range")
    parser.add_argument("-c", "--cargo-concentration", type=float, help="the concentration of cargo in "\
                        "the cell, in Molars", default=50e-6)
    parser.add_argument("-pkl", "--pickle-file", type=str, \
                        help="Filename of output pickle file. If supplied, the heatmap will be"\
                             " pickled and saved to this file.")
    parser.add_argument("-o", "--output", type=str, help="Output filename prefix. The final output"\
                        " filenames will be {output}_p_{passive}.png, for each passive in "\
                        "--passive-range")
    args = parser.parse_args()

    if args.output is None:
        args.output = f"Heatmaps_for_c_{args.cargo_concentration}_"\
                      f"n_{args.npc_traverse_range[0]}_{args.npc_traverse_range[1]}_"\
                      f"K_{args.k_on_range[0]}_{args.k_on_range[1]}"
    print(args)
    main(args.output, args.passive_range, args.npc_traverse_range, args.k_on_range, \
         args.nx, args.ny, args.n_passive, args.cargo_concentration, args.pickle_file)
