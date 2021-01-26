# ---
# jupyter:
#   jupytext:
#     formats: py:hydrogen
#     text_representation:
#       extension: .py
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.5.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import multiprocessing
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import transport_simulation 

# %%


def get_tau_passive_diffusion(nmol_per_sec_per_M, volume_L):
    ''' 
    This utility function computes tau for one-sided passive diffusion
    specified in number of molecules per second per M, given
    the volume from which the passive diffusion leaves
    TODO: is this needed anywhere?
    '''
    global N_A
    gamma= nmol_per_sec_per_M/N_A/volume_L
    tau= 1.0/gamma
    return tau


def get_new_transport_simulation(**kwargs):
    ts= transport_simulation.TransportSimulation(**kwargs) 
    ts.set_time_step(1e-3)
    return ts

def mp_do_simulation(param_range, i, j, 
                     equilibration_time_sec,
                    transport_simulation_generator):

    nskip_statistics= 100
    my_ts= transport_simulation_generator()
    cur_params= { param_range['tag_x']: param_range['range_x'][i],
                  param_range['tag_y']: param_range['range_y'][j]}
    my_ts.set_params(**cur_params)
    stats = my_ts.simulate(equilibration_time_sec,
                           nskip_statistics= 100)  
    if (np.random.rand()<0.1):
        print(f"Finished some ten jobs (at i={i} j={j})")
    return {"i":i, "j":j, "stats":stats}

def mp_handle_stats(stats_grids, mydicts):
    for mydict in mydicts:
        i= mydict["i"]
        j= mydict["j"]
        stats= mydict["stats"]
#        print("My Handle - ", i, j)
        for nmol_type in stats.keys():
            if nmol_type=="time_sec":
                continue
            try: 
                pass
                stats_grids[nmol_type][j, i]= stats[nmol_type][-1] # matrix is indexed row first 
            except KeyError:
                print(f"my_handle_stats - Key {nmol_type} not found")
            except IndexError as e:
                print(e)

def mp_handle_error(error):
    print("Error", error)
    
def map_param_grid_parallel(param_range, 
                   equilibration_time_sec= 150.0, 
                   n_processors= 5,
                   transport_simulation_generator= get_new_transport_simulation
                  ):
    '''
    Run in parallel to compute simulation results over a range of parameters
    
    :param param_range: a dictionary with 'tag_x'/'tag_y' keys corresponding
                 to x/y-axis param names, resp., and 'range_x'/'range_y' keys
                 for numpy array for corresponding param ranges
    :param equilibration_time_sec: equilibration time per condition
    :param n_processors: number of processors on which to run in parallel
    
    :return a dictionary from simulation properties to 2D arrays with the
            their values at the end of simulations for each parameter
            combination
    '''
    VERBOSE= True
    nx= len(param_range['range_x'])
    ny= len(param_range['range_y'])
    ts= transport_simulation_generator() # a sample ts to be returned
    if VERBOSE:
        for param in [param_range['tag_x'], param_range['tag_y']]:
            print("Param {:} default value is {:}".format(param, getattr(ts, param)))
    stats_grids = {}
    for nmol_type in ts.nmol.keys(): # TODO: add get_nmols()
        stats_grids[nmol_type]= np.ndarray((ny, nx)) # Row major - y coordinate goes first
    jobs_params= []
    for j in range(ny):
        for i in range(nx):
            jobs_params.append((param_range.copy(),
                                i,
                                j,
                                equilibration_time_sec,
                               transport_simulation_generator))
    print("njobs={}".format(len(jobs_params)))
    callback_function = \
        lambda mydict: mp_handle_stats(stats_grids, mydict)
    pool= multiprocessing.Pool(processes= n_processors)
    results= pool.starmap_async(mp_do_simulation,
                               jobs_params,
                               callback= callback_function,
                               error_callback=mp_handle_error)
    results.wait()
    pool.close()
    pool.join()
    return stats_grids, ts



###### VISUALIZATION OF RESULTS ######

def plot_param_grid(param_range, 
                    Z, 
                    Z_label= None, 
                    **contourf_kwargs):
    x_meshgrid, y_meshgrid = np.meshgrid(param_range["range_x"],
                                         param_range["range_y"])
    plt.contourf( x_meshgrid, 
                  y_meshgrid, 
                  Z,
                **contourf_kwargs)
    ax= plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(param_range['pretty_x'])
    ax.set_ylabel(param_range['pretty_y'])
    xlim= ax.get_xlim()
    ylim= ax.get_ylim()
    if xlim[1]>xlim[0]:
        ax.set_xlim(xlim[1], xlim[0])
    cb = plt.colorbar(label= Z_label)
    ticks = cb.get_ticks()
    cb.set_ticks(ticks)

def get_N_to_C_ratios(stats_grids, v_N_L, v_C_L):
    ''' return N/C ratios from stats_grids computed in the previous cell'''
    nNs= stats_grids["complexL_N"]+stats_grids["freeL_N"]+stats_grids["complexU_N"]+stats_grids["freeU_N"] 
    nCs= stats_grids["complexL_C"]+stats_grids["freeL_C"]+stats_grids["complexU_C"]+stats_grids["freeU_C"]
    ratios= (nNs/v_N_L) / (nCs/v_C_L)
    return ratios

def plot_NC_ratios(param_range, stats_grids, ts, 
                   vmin= 1.0, vmax= 4.0, 
                   locator= None, 
                  levels= None):
    NC_ratios= get_N_to_C_ratios(stats_grids, 
                              v_N_L= ts.get_v_N_L(), 
                              v_C_L= ts.get_v_C_L())
    print(vmin, vmax)
    plot_param_grid(param_range, 
                    NC_ratios,
                    Z_label= "N/C ratio",
                    vmin= vmin,
                    vmax= vmax, 
                    locator= locator, 
                    levels= levels,
                    extend='both')

def plot_bound_fraction(param_range, stats_grids, compartment, ax= None):
    assert(compartment in ['N', 'C'])
    tag_complex= f"complexL_{compartment}"
    tag_free= f"freeL_{compartment}"
    complexL_fraction= stats_grids[tag_complex]/(stats_grids[tag_complex]+stats_grids[tag_free])
    ##### Contourf
    vmax=1.0
    if ax != None:
        plt.sca(ax)
    plot_param_grid(param_range, 
                    complexL_fraction,
                    Z_label= f"bound fraction ({compartment})",
                    vmin=0.0, 
                     vmax=vmax, 
                     levels=np.linspace(0,vmax,11), 
                     extend='both')

def get_import_export_ratios(stats_grids):
    ''' return import/export ratios from stats_grids computed in the previous cell'''
    import_rate= stats_grids["nuclear_importL_per_sec"]+stats_grids["nuclear_importU_per_sec"]
    export_rate= stats_grids["nuclear_exportL_per_sec"]+stats_grids["nuclear_exportU_per_sec"]
    ratios= import_rate/export_rate
    return ratios
    
def plot_import_export(param_range, stats_grids, axes= None,
                      **contourf_kwargs):
    if axes==None:
        fig, axes= subplot(2, 1)
    for tag, ax in zip(['import', 'export'], axes):
        plt.sca(ax)
        plot_param_grid(param_range, 
                        stats_grids[f'nuclear_{tag}L_per_sec'],
                        Z_label= tag + r' [$sec^{-1}$]',
                       **contourf_kwargs)

# %%
