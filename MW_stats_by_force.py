
import multiprocessing
import transport_simulation 
import sys
from transport_simulation import TransportSimulation

__all__ = ["get_MW_stats_list_by_force"]

s = 20.0
no_force_coefficient = 3.0
force_coefficient = 10.0

def do_simulate(ts, simulation_time_sec):
    return ts.simulate(simulation_time_sec)

def get_ts_with_parameters(MW= 27, 
                      NLS_strength= 0, 
                      is_force= False, 
                      **kwargs):
    if is_force:
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
                54:0.00026308, 
                67:0.00272423 }
    return effects[MW]

def get_fraction_complex_NPC_traverse_per_sec(MW, is_force):
    rate= { 27: [s*no_force_coefficient,  s*force_coefficient],
            41: [s*no_force_coefficient,  s*force_coefficient],
            54: [s*no_force_coefficient, s*force_coefficient],
            67: [s*no_force_coefficient,  s*force_coefficient] }
    i_force= 1 if is_force else 0
    return rate[MW][i_force]

def get_MW_stats_list_by_force(MW, simulation_time_sec,
                              n_processors=None, nls_range=(0,9)):
    assert(MW in [27,41, 54, 67])
    if n_processors is None:
        n_processors= multiprocessing.cpu_count()
        
    stats_list_by_force= {}
    TSs_by_force= {}
    for is_force in [False, True]:
        TS_tuples= []
        for i_NLS in range(*nls_range):
            ts = get_ts_with_parameters(MW= MW,
                                    NLS_strength=i_NLS,
                                  is_force= is_force)
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
        s = float(sys.argv[3])
        no_force_coefficient = float(sys.argv[3])
        force_coefficient = float(sys.argv[4])
    else:
        s = 20.0
        no_force_coefficient = 3.0
        force_coefficient = 10.0
    
    result = get_MW_stats_list_by_force(MW, simulation_time_sec)
    import pickle
    with open(f"MW_stats_list_{MW}_{simulation_time_sec}_{s*no_force_coefficient}_{s*force_coefficient}.pkl", 'wb') as f:
        pickle.dump(result, f)




