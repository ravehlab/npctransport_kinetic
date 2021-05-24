
import sys
import multiprocessing
import transport_simulation 
from transport_simulation import TransportSimulation
from make_plots import make_plot

__all__ = ["get_MW_stats_list_by_force"]

no_force = 30.0
force = 200.0

def do_simulate(ts, simulation_time_sec):
    return ts.simulate(simulation_time_sec)

def get_ts_with_parameters(MW= 27, 
                           NLS_strength= 0, 
                           is_force= False, 
                           is_change_cell_volume= True,
                           **kwargs):
    if is_force and is_change_cell_volume:
        v_N_L=762e-15
        v_C_L=4768e-15
    else:
        v_N_L=627e-15
        v_C_L=2194e-15    
    ts= transport_simulation.TransportSimulation(v_N_L= v_N_L,
                                                v_C_L= v_C_L)
    ts.set_time_step(0.1e-3)
    ts.set_NPC_dock_sites(n_NPCs= 2000, 
                        n_dock_sites_per_NPC= 500)
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
    #
    ts.set_passive_nuclear_molar_rate_per_sec(
        get_passive_nuclear_molar_rate_per_sec(MW, is_force))
    ts.set_params(rate_free_to_complex_per_sec= 
                  get_free_to_complex_rate(NLS_strength))
    ts.set_params(fraction_complex_NPC_traverse_per_sec=
                  get_fraction_complex_NPC_traverse_per_sec(MW, is_force))
    #
    ts.set_params(**kwargs) # override defaults
    return ts

def get_free_to_complex_rate(NLS_strength):
    rates = [0.0,
             0.001,
             0.00316,
             0.01,
             0.02, #2.11
             0.045, #2.11
             0.1,  #16.4
             0.2,
             0.45,
             1.0,
             2.0,
             4.5
            ]
    return rates[NLS_strength]

def get_passive_nuclear_molar_rate_per_sec(MW, is_force): # TODO: verify it corresponds to multiplyng by concentration rather than nmolecules
    #TODO: generalize this - either from the literature or regression
    base_rates={ 27:0.0805618, 
                41:0.06022355, 
                54:0.03301662, 
                67:0.0287649 }
    rate= base_rates[MW]
    if is_force:
        rate += get_force_effect_on_diffusion(MW)
    return rate

def get_force_effect_on_diffusion(MW):
    """
    The effect of force on passive diffusion as measured by experiment
    """
    effects = {27:0.08214946, 
                41:0.03027974, 
                54:0.01, # 54:0.00026308, 
                67:0.01 } #67:0.00272423 }
    return effects[MW]

def get_fraction_complex_NPC_traverse_per_sec(MW, is_force):
    rate_row = [no_force, force]
    rate= { mw : rate_row for mw in [27, 41, 54, 67] }
    i_force= 1 if is_force else 0
    return rate[MW][i_force]

def get_MW_stats_list_by_force(MW, simulation_time_sec, n_processors=None, \
                               is_change_cell_volume=True, nls_range=(0,9)):
    assert(MW in [27,41, 54, 67])
    if n_processors is None:
        n_processors= multiprocessing.cpu_count()
        
    stats_list_by_force= {}
    TSs_by_force= {}
    for is_force in [False, True]:
        TS_tuples= []
        for i_NLS in range(*nls_range):
            ts = get_ts_with_parameters(MW=MW,
                                        NLS_strength=i_NLS,
                                        is_force=is_force,
                                        is_change_cell_volume=is_change_cell_volume)
            TS_tuples.append((ts, simulation_time_sec))
        pool= multiprocessing.Pool(processes= n_processors)
        stats_list_by_force[is_force]= pool.starmap(do_simulate,
                                                    TS_tuples)
        TSs_by_force[is_force]= [x[0] for x in TS_tuples]
        print(f"Is force {is_force} i_NLS {i_NLS}: OK")
    return (stats_list_by_force, TSs_by_force)

if __name__ == "__main__":
    MW = int(sys.argv[1])
    simulation_time_sec = float(sys.argv[2])
    
    if len(sys.argv) > 3:
        no_force = float(sys.argv[3])
        force = float(sys.argv[4])
    else:
        no_force = 30.0
        force = 200.0
    filename = f"MW_stats_list_{MW}_{simulation_time_sec}_{no_force}_{force}"
    
    # result: is_force -> [stats_dictionary for NLS in free_to_complex_rates]
    result = get_MW_stats_list_by_force(MW, simulation_time_sec)
    make_plot(result, f"{filename}.png")
    print("Figure saved as {filename}.png")

    import pickle

    # the keys of all the values we need for the final graphs:
    keys_for_graphs = ['nuclear_importL_per_sec', 'nuclear_importU_per_sec', 
                       'nuclear_exportL_per_sec', 'nuclear_exportU_per_sec']
    for label in ['L', 'U']:
        for side in ['N', 'C']:
            # nmol in NPC
            for source in ['import', 'export']:
                keys_for_graphs.append('complex{}_NPC_{}_{}'.format(label, side, source))
            # nmol in nucleus and cytoplasm
            for state in ['free','complex']:
                keys_for_graphs.append('{}{}_{}'.format(state, label, side))

    final_result = dict()
    for key in keys_for_graphs:
        final_result[key] = {}
        # TODO: double check this
        for is_force in [False, True]:
            for i_NLS, stats in enumerate(result[is_force]):
                final_result[key][(i_NLS, is_force)] = stats[key][-1]

    with open(f"final_{MW}_{simulation_time_sec}_{no_force}_{force}.pkl", 'wb') as f:
        pickle.dump(final_result, f)
        print(f"Saved final results")

    with open(f"{filename}.pkl", 'wb') as f:
        pickle.dump(result, f)
        print(f"Saved as {filename}.pkl")



