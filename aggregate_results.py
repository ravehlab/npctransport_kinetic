import sys
import pickle
import numpy as np
import transport_simulation



def get_compartment_concentration_stats(stats, compartment, v_N_L=627e-15, v_C_L=2194e-15):
    assert(compartment in ['N','C'])
    nmol_stats= get_compartment_nmol_stats(stats, compartment)
    is_nuclear= (compartment=='N')
    volume_L= (v_N_L if is_nuclear else v_C_L)
    return (nmol_stats/transport_simulation.N_A)/volume_L

def get_compartment_nmol_stats(stats, compartment):
    nmol_stats= np.zeros((2, nNLS))
    for state in ['free','complex']:
        for label in ['L', 'U']:
            tag= '{}{}_{}'.format(state, 
                                  label, 
                                  compartment)
            nmol_stats[0] += [stats[tag][i, False] for i in range(nNLS)]
            nmol_stats[1] += [stats[tag][i, True] for i in range(nNLS)]
    return nmol_stats

def get_N2C_ratios(stats):
    nmol_N = get_compartment_concentration_stats(stats, "N")
    nmol_C = get_compartment_concentration_stats(stats, "C")
    ratios = nmol_N/nmol_C
    return ratios[0], ratios[1]

def get_export_rates(stats):
    """get_export_rates.

    Parameters
    ----------
    stats :
        stats dictionary resulting from unpacking the pickled final_MW_[...].pkl file from 
        MW_stats_by_force

    Returns
    -------
    exp_rates_no_force : list
        exp_rates_no_force[i] is the export rate with i_NLS and no force
    exp_rates_force : list
        exp_rates_force[i] is the export rate with i_NLS and force  

    """
    exp_rates_no_force = [stats['nuclear_exportL_per_sec'][i, False] + \
                          stats['nuclear_exportU_per_sec'][i, False] for i in range(nNLS)]

    exp_rates_force = [stats['nuclear_exportL_per_sec'][i, True] + \
                       stats['nuclear_exportU_per_sec'][i, True] for i in range(nNLS)]

    return exp_rates_no_force, exp_rates_force

def get_import_rates(stats):
    """get_import_rates.

    Parameters
    ----------
    stats :
        stats dictionary resulting from unpacking the pickled final_MW_[...].pkl file from 
        MW_stats_by_force

    Returns
    -------
    imp_rates_no_force : list
        imp_rates_no_force[i] is the import rate with i_NLS and no force
    imp_rates_force : list
        imp_rates_force[i] is the import rate with i_NLS and force  

    """
    imp_rates_no_force = [stats['nuclear_importL_per_sec'][(i, False)] + \
                          stats['nuclear_importU_per_sec'][(i, False)] for i in range(nNLS)]

    imp_rates_force = [stats['nuclear_importL_per_sec'][(i, True)] + \
                       stats['nuclear_importU_per_sec'][(i, True)] for i in range(nNLS)]

    return imp_rates_no_force, imp_rates_force

if __name__ == "__main__":
    sec = float(sys.argv[1])
    no_force = float(sys.argv[2])
    force = float(sys.argv[3])
    results = {}
    for MW in [27, 34, 41, 47, 54, 67]:
        with open(f"final_{MW}_{sec}_{no_force}_{force}.pkl", 'rb') as f:
            results[MW] = pickle.load(f)

    nNLS = len(results[27]['freeL_C'])//2

    final_values = {}

    final_values['import'] = {}
    final_values['export'] = {}
    final_values['N2C'] = {}
    for MW in results.keys():
        final_values['import'][MW] = get_import_rates(results[MW])
        final_values['export'][MW] = get_export_rates(results[MW])
        final_values['N2C'][MW] = get_N2C_ratios(results[MW])

    with open(f"stats_for_graph_{no_force}_{force}.pkl", 'wb') as f:
        pickle.dump(final_values, f)
