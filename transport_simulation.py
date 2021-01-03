
import numpy as np
import re

N_A= 6.022e+23 # Avogadro's number 

# TODO: test
def fL_to_L(v_fL):
    ''' convert femtoliters to liters '''
    return v_fL*1E-15

def register_update_functions(cls):
    """
    A class decorator that allows registering methods. 
    Methods registered will be added to a class variable list called '_update_funcs'
    To register a function, use the decorator 'register_update()' defined below.
    """
    cls._update_funcs = []
    for methodname in dir(cls):
        method = getattr(cls, methodname)
        if hasattr(method, '_update_func'):
            if method._update_func:
                cls._update_funcs.append(method)
    return cls

def register_update(active=True):
    """
    A function which created a decorator that adds a function a bool attribute called '_update_func'
    with the value 'active'
    """
    def wrapper(func):
        func._update_func = active
        return func
    return wrapper

def move_nmol(T_list, src, dst, nmol):
    T_list.append((src, dst, nmol))

@register_update_functions
class TransportSimulation():

    ###########################################
    # getter/setter functions / utility
    ###########################################

    def set_nmol(self, species, value):
        self.nmol[species]= value

    def get_nmol(self, species):
        return self.nmol[species]

    def get_compartment(self, species):
        m= re.search("_([a-zA-Z]*)$", species)
        if m is None:
            raise ValueError(f"Can't parse compartment from species {species}")
        return m.group(1)

    def get_compartment_volume_L(self, species):
        compt= self.get_compartment(species)
        if compt=="N":
            return self.v_N_L
        elif compt=="C":
            return self.v_C_L
        elif compt=="cell":
            return self.v_N_L + self.v_C_L
        else:
            raise ValueError(f"Only nucleus/cytoplasm/cell has a volume (species {species} compartment {compt})") 

    def set_concentration_M(self, species, c_M):
        ''' 
        Sets the concentration of specified speciecs to c_M (in M units) 
        @raise ValueError if compartment has no volume
        '''
        global N_A
        v_L= self.get_compartment_volume_L(species)
        self.nmol[species]= c_M * v_L * N_A

    def get_concentration_M(self, species):
        ''' 
        returns the concentrations of specied species in its compartment 
        @raise ValueError if compartment has no volume
        '''
        global N_A
        v_L= self.get_compartment_volume_L(species)
        return self.nmol[species]/(N_A*v_L)

    def set_RAN_distribution(self,
                             Ran_cell_M,
                             parts_GTP_N,
                             parts_GTP_C,
                             parts_GDP_N,
                             parts_GDP_C):
        '''
        Sets the RAN distribution among compartments based on relative parts specified, s.t. total concentration is constant

        @param Ran_cell_M - total Ran concentration in cell
        @param GTP_N, GTP_C, GDP_N, GDP_C - relative quantities of all Ran species in nucleus (N) and cytoplasm (N)
        '''
        global N_A
        RAN_distribution = np.array([parts_GTP_N, parts_GTP_C,
                                     parts_GDP_N, parts_GDP_C]) 
        RAN_distribution = RAN_distribution/np.sum(RAN_distribution) # normalize to 1
        nmol_Ran_cell= Ran_cell_M * (self.v_N_L + self.v_C_L) * N_A 
        self.nmol["GTP_N"] = nmol_Ran_cell * RAN_distribution[0]
        self.nmol["GTP_C"] = nmol_Ran_cell * RAN_distribution[1]
        self.nmol["GDP_N"] = nmol_Ran_cell * RAN_distribution[2]
        self.nmol["GDP_C"] = nmol_Ran_cell * RAN_distribution[3]

    def set_v_N_L(self, v_L):
        self.v_N_L= v_L

    def set_v_C_L(self, v_L):
        self.v_C_L= v_L

    def get_v_N_L(self):
        return self.v_N_L

    def get_v_C_L(self):
        return self.v_C_L

    def get_v_cell_L(self):
        return self.v_N_L + self.v_C_L

    def set_time_step(self, dt_sec):
        ''' set time step in seconds '''
        self.dt_sec= dt_sec

    def get_time_step(self, dt_sec):
        ''' get time step in seconds '''
        self.dt_sec= dt_sec  

    def reset_simulation_time(self):
        self.sim_time_sec= 0.0


    ###################
    # Consturctor (and init functions)
    ###################
    def set_passive_diffusion_rate_per_sec(self, rate):
        '''
        Sets the parameter max_passive_diffusion_rate_nmol_per_sec_per_M such that the percent of 
        molecules that passively diffuse per second is the argument 'rate'
        '''
        self.max_passive_diffusion_rate_nmol_per_sec_per_M = rate*N_A
    def set_params(self, **kwargs):
        for param, value in kwargs.items():
            assert hasattr(self, param)
            setattr(self, param, value)

    def _init_simulation_parameters(self, **kwargs):
        # TODO: add all simulation parameters here with proper units
        self.dt_sec = 1e-4 # simulation time step  
        # NPC dock capacity:
        n_NPCs= 200 # (maximal estimate from Timney et al. 2016 paper)
        n_dock_sites_per_NPC= 500 #  dock sites for cargo-importin complexes per NPC, rule of thumb estimate  # TODO: this may depend on molecule size
        self.NPC_dock_sites = n_NPCs * n_dock_sites_per_NPC # total capacity for cargo-importin complexes in entire NPC, in number of molecules
        # Rates:  # TODO: change nmol to nmolec - to prevent confusion between moles and molecules
        self.rate_complex_to_NPC_per_free_site_per_sec_per_M = 50000.0e+6/self.NPC_dock_sites # the fraction of cargo-importin complexes that will dock to avaialble NPC dock sites per second (from either cytoplasm or nucleus)
        self.fraction_complex_NPC_to_free_N_per_M_GTP_per_sec = 0.005e+6
        self.fraction_complex_N_to_free_N_per_M_GTP_per_sec = 0.005e+6
        self.fraction_complex_NPC_to_complex_N_C_per_sec= 1.0
        self.rate_GDP_N_to_GTP_N_per_sec= 200.0
        self.rate_GTP_N_to_GDP_N_per_sec= 0.2
        self.rate_GTP_C_to_GDP_C_per_sec= 500.0
        self.rate_GTP_N_to_GTP_C_per_sec = 0.15
        self.rate_GDP_C_to_GDP_N_per_sec = 0.2
        self.rate_GDP_N_to_GDP_C_per_sec = 0.2
        self.rate_complex_C_to_free_C_per_sec = 0.05       
        self.rate_free_to_complex_per_sec = 0.10 # assuming importins are not rate limiting in either cytoplasm aor nucleus and have identical concentration
        self.passive_competition_weight= 0.0 # a number between 0.0 and 1.0 quantifying the weight of competition # TODO: this could be a flag
        self.max_passive_diffusion_rate_nmol_per_sec_per_M= 20000 # as the name suggests, without accounting for competition effects # TODO: in future, a single number for both import and export that is independent of C/N volumes, # of NPCs etc
        self.bleach_volume_L_per_sec= 1.0e-15 # cytoplasmic cargo volume being bleached per second
        self.bleach_start_time_sec= np.inf # no bleaching by default
        self.set_params(**kwargs)


    def __init__(self, **kwargs):
        ''' Set initial state of the simulation '''
        self._init_simulation_parameters(**kwargs)
        self.sim_time_sec= 0.0
        self.nmol= {} # number of molecules of various species
        # Cell geometry:
        self.v_C_L= 10e-15 # Cytoplsmic volume in L
        self.v_N_L= 3e-15 # Nuclear volume in L
        # NPC:
        self.nmol["complexL_NPC"]= 1e0 # number of cargo-importin complexes docked to the NPC (labeled)
        self.nmol["complexU_NPC"]= 0 # (unlabeled)
        # NPC directionality:
        self.NPC_N2C_ratio = 0.5 # ratio of cargo in the NPC that docked from the nucleus
        # Cytoplasm:
        self.set_concentration_M("cargo_C", 50e-6)  # Nuclear concentration of labeled cargo in M
        self.nmol["complexL_C"]=  self.get_nmol("cargo_C")*0.25 # number of cargo-importin complexes in cytoplasm (labeled)
        self.nmol["freeL_C"]= self.nmol["cargo_C"] - self.nmol["complexL_C"] # number of free cargo molecules in cytoplasm (labeled)
        self.nmol["complexU_C"]= 0 # (unlabeled)
        self.nmol["freeU_C"]= 0 # (unlabeled)
        del self.nmol["cargo_C"]
        # Nucleus:
        self.set_concentration_M("cargo_N", 0e-5)  # Nuclear concentration of labeled cargo in M
        self.nmol["complexL_N"] = 0 # number of cargo-importin complexes in nucleus (labeled)
        self.nmol["freeL_N"]= self.nmol["cargo_N"] - self.nmol["complexL_N"] # number of free cargo molecules in nucleus (labeled)
        self.nmol["complexU_N"]= 0 # (unlabeled)
        self.nmol["freeU_N"]= 0 # (unlabeled)
        del self.nmol["cargo_N"] 
        # import export per dt_sec
        self.nmol["import"] = 0
        self.nmol["export"] = 0
        # Ran in all:
        self.set_RAN_distribution(Ran_cell_M= 2e-5, # total physiological concentration of Ran # TODO: check in the literature 
                                  parts_GTP_N=1000,
                                  parts_GTP_C=1,
                                  parts_GDP_N=1,
                                  parts_GDP_C=1000)

    ##########################
    # Transitions calculators:
    ########################

    @register_update()
    def get_nmol_complex_NPC_to_free_N(self, T_list):
        """
        Number of labeled cargo molecules released from the NPC to the nucleus over a self.dt_sec time step
        (Note: it is assumed each undocking leads to export of a single RanGTP molecule)

        Return: dictionary with number of molecules to add/subtract from each species
        """
        #return float(int(np.power(nmol_GTP_N/max_RAN, 5)*nmol_NPC))
        f= self.fraction_complex_NPC_to_free_N_per_M_GTP_per_sec  \
            * self.get_concentration_M("GTP_N") \
            * self.dt_sec
        nL= f * self.nmol["complexL_NPC"] 
        nU= f * self.nmol["complexU_NPC"] 
        n= nL+nU

        assert n <= self.nmol["GTP_N"] and nL <= self.nmol["complexL_NPC"] and nU <= self.nmol["complexU_NPC"]            

        move_nmol(T_list,
                   src="complexL_NPC",\
                   dst="freeL_N",\
                   nmol=nL)
        move_nmol(T_list,
                   src="complexU_NPC",\
                   dst="freeU_N",\
                   nmol=nU)
        move_nmol(T_list,
                   src="GTP_N",\
                   dst="GTP_C",\
                   nmol=n)


    @register_update()
    def get_nmol_complex_N_to_free_N(self, T_list):
        """
        Number of labeled cargo molecules that disassemble in the nucleus over a self.dt_sec time step
        Note: it is assumed each undocking leads to export of a single RanGTP molecule instantaneously

        Return: dictionary with number of molecules to add/subtract from each species
        """
        #return float(int(np.power(nmol_GTP_N/max_RAN, 5)*nmol_NPC))
        c_GTP_N_M= self.get_concentration_M("GTP_N")
        f= self.fraction_complex_N_to_free_N_per_M_GTP_per_sec \
            * self.get_concentration_M("GTP_N") \
            * self.dt_sec
        nL= f * self.nmol["complexL_N"] 
        nU= f * self.nmol["complexU_N"] 
        n= nL+nU                      
        #     print("n {} GTP_N {} complex_N {}".format(n, self.nmol["GTP_N"], self.nmol["complex_N"]))
        assert n <= self.nmol["GTP_N"] and nL <= self.nmol["complexL_N"] and nU <= self.nmol["complexU_N"]            
        move_nmol(T_list,
                   src="complexL_N",\
                   dst="freeL_N",\
                   nmol=nL)
        move_nmol(T_list,
                   src="complexU_N",\
                   dst="freeU_N",\
                   nmol=nU)
        move_nmol(T_list,
                   src="GTP_N",\
                   dst="GTP_C",\
                   nmol=n)

    @register_update()
    def get_nmol_GDP_N_to_GTP_N(self, T_list):
        """
        Number of GDP molecules in the nucleus converted to GTP

        Return: dictionary with number of molecules to add/subtract from each species
        """
        n1= self.rate_GDP_N_to_GTP_N_per_sec \
            * self.nmol["GDP_N"] \
            * self.dt_sec
        n2= self.rate_GTP_N_to_GDP_N_per_sec \
            * self.nmol["GTP_N"] \
            * self.dt_sec
        n= n1-n2
        move_nmol(T_list,
                   src="GDP_N",\
                   dst="GTP_N",\
                   nmol=n)

    @register_update()
    def get_nmol_GTP_C_to_GDP_C(self, T_list):
        """
        Number of GTP molecules in the cytoplasm converted to GDP

        Return: dictionary with number of molecules to add/subtract from each species
        """
        n= self.rate_GTP_C_to_GDP_C_per_sec \
            * self.nmol["GTP_C"] \
            * self.dt_sec
        move_nmol(T_list,
                   src="GTP_C",\
                   dst="GDP_C",\
                   nmol=n)

    @register_update()
    def get_nmol_GTP_N_to_GTP_C(self, T_list):
        """
        Number of GTP molecules exported from the nucleus

        Return: dictionary with number of molecules to add/subtract from each species
        """

        n= self.rate_GTP_N_to_GTP_C_per_sec \
            * self.nmol["GTP_N"] \
            * self.dt_sec
        move_nmol(T_list,
                   src="GTP_N",\
                   dst="GTP_C",\
                   nmol=n)

    @register_update()
    def get_nmol_GDP_C_to_GDP_N(self, T_list):
        """
        Number of GDP molecules imported to the nucleus

        Return: dictionary with number of molecules to add/subtract from each species
        """
        n1= self.rate_GDP_C_to_GDP_N_per_sec \
            * self.nmol["GDP_C"] \
            * self.dt_sec
        n2= self.rate_GDP_N_to_GDP_C_per_sec \
            * self.nmol["GDP_N"] \
            * self.dt_sec
        n= n1-n2
        move_nmol(T_list,
                   src="GDP_C",\
                   dst="GDP_N",\
                   nmol=n)

    @register_update()
    def get_nmol_complex_C_to_free_C(self, T_list):
        """
        The number of cargo-importin complexes that unbind importin over time step dt_sec

        Return: dictionary with number of molecules to add/subtract from each species
        """
        f= self.rate_complex_C_to_free_C_per_sec \
            * self.dt_sec
        nL= f * self.nmol["complexL_C"]
        nU= f * self.nmol["complexU_C"]
        assert(nL <= self.nmol["complexL_C"])
        assert(nU <= self.nmol["complexU_C"])
        move_nmol(T_list,
                   src="complexL_C",\
                   dst="freeL_C",\
                   nmol=nL)
        move_nmol(T_list,
                   src="complexU_C",\
                   dst="freeU_C",\
                   nmol=nU)

    @register_update()
    def get_nmol_free_C_to_complex_C(self, T_list): # assume importin is not rate limiting
        """
        The number of the labeled molecules that bind to importin over time step dt_sec
        in the cytoplasm

        Return: dictionary with number of molecules to add/subtract from each species
        """
        f= self.rate_free_to_complex_per_sec \
            * self.dt_sec
        nL= f * self.nmol["freeL_C"]
        nU= f * self.nmol["freeU_C"]
        assert(nL <= self.nmol["freeL_C"])
        assert(nU <= self.nmol["freeU_C"])
        move_nmol(T_list,
                   src="freeL_C",\
                   dst="complexL_C",\
                   nmol=nL)
        move_nmol(T_list,
                   src="freeU_C",\
                   dst="complexU_C",\
                   nmol=nU)

    @register_update()
    def get_nmol_free_N_to_complex_N(self, T_list): # assume importin is not rate limiting
        """
        The number of the labeled molecules that bind to importin over time step dt_sec
        in the nucleus

        Return: dictionary with number of molecules to add/subtract from each species
        """
        f= self.rate_free_to_complex_per_sec \
            * self.dt_sec
        nL= f * self.nmol["freeL_N"]
        nU= f * self.nmol["freeU_N"]
        assert(nL <= self.nmol["freeL_N"])
        assert(nU <= self.nmol["freeU_N"])
        move_nmol(T_list,
                   src="freeL_N",\
                   dst="complexL_N",\
                   nmol=nL)
        move_nmol(T_list,
                   src="freeU_N",\
                   dst="complexU_N",\
                   nmol=nU)


    @register_update()
    def get_free_N_to_free_C(self, T_list): # passive
        """
        Computes the net number of unbound molecules in the nucleus that passively export 
        to the cytoplasm per second (net = export - import)

        Return: dictionary with number of molecules to add/subtract from each species

        # COMMENT: a proper treatment of this would depend on ratio between nuclear 
        # and cytoplasmic volumes, number of NPCs etc - here we ignore this for now
        # - we can change it in future based on theoretical equations of passive diffusion
        """
        # Comment: competition is assumed to have zero effect at this time
        fraction_bound_dock_sites_NPC= (self.nmol["complexL_NPC"] + self.nmol["complexU_NPC"]) / self.NPC_dock_sites
        competition_multiplier= 1.0 - self.passive_competition_weight * fraction_bound_dock_sites_NPC
        f= self.max_passive_diffusion_rate_nmol_per_sec_per_M \
            * competition_multiplier \
            * self.dt_sec

        nL_export = f * self.get_concentration_M("freeL_N")  
        nL_import = f * self.get_concentration_M("freeL_C")
        nU_export = f * self.get_concentration_M("freeU_N")  
        nU_import = f * self.get_concentration_M("freeU_C")
        
        move_nmol(T_list,
                   src="freeL_N",\
                   dst="freeL_C",\
                   nmol=nL_export)
        move_nmol(T_list,
                   src="freeU_N",\
                   dst="freeU_C",\
                   nmol=nU_export)
        move_nmol(T_list,
                   src="freeL_C",\
                   dst="freeL_C",\
                   nmol=nL_import)
        move_nmol(T_list,
                   src="freeU_C",\
                   dst="freeU_N",\
                   nmol=nU_import)

    @register_update()
    def get_nmol_complex_N_C_to_complex_NPC(self, T_list):
        """
        Computes the number of molecules that bind to the NPC from the nucleus
        over dt_sec time step (These will all be bound to importin)

        Return: dictionary with number of molecules to add/subtract from each species
        """
        nmol_free_sites_NPC = (self.NPC_dock_sites - self.nmol["complexL_NPC"] - self.nmol["complexU_NPC"])
        f = nmol_free_sites_NPC \
            * self.rate_complex_to_NPC_per_free_site_per_sec_per_M \
            * self.dt_sec
        # TODO: debug - something is weird here (BR Dec 11,2020)
        cL_N_M= self.get_concentration_M("complexL_N")
        cL_C_M= self.get_concentration_M("complexL_C")
        cU_N_M= self.get_concentration_M("complexU_N")
        cU_C_M= self.get_concentration_M("complexU_C")
        nL_N= f * cL_N_M
        nL_C= f * cL_C_M 
        nU_N= f * cU_N_M 
        nU_C= f * cU_C_M
        assert_coeff= 2.0
        assert1_almost= (nL_N+nL_C+nU_N+nU_C <= assert_coeff*nmol_free_sites_NPC) 
        assert2_almost= (nL_N+nU_N <= assert_coeff*(self.nmol["complexL_N"]+self.nmol["complexU_N"])) 
        assert3_almost= (nL_C+nU_C <= assert_coeff*(self.nmol["complexL_C"]+self.nmol["complexU_C"]))
        if(not (assert1_almost and assert2_almost and assert3_almost)):       
            assert1= (nL_N+nL_C+nU_N+nU_C <= nmol_free_sites_NPC) 
            assert2= (nL_N+nU_N <= self.nmol["complexL_N"]+self.nmol["complexU_N"]) 
            assert3= (nL_C+nU_C <= self.nmol["complexL_C"]+self.nmol["complexU_C"])
            print(self.nmol)
            print(f"f {f} dLabeled: N {nL_N} C {nL_C}, dUnlabeled: N {nU_N} C {nU_C}")
            assert(assert1)
            assert(assert2)
            assert(assert3)
            assert(assert1 and assert2 and assert3)
        move_nmol(T_list,
                   src="complexL_N",\
                   dst="complexL_NPC",\
                   nmol=nL_N)
        move_nmol(T_list,
                   src="complexL_C",\
                   dst="complexL_NPC",\
                   nmol=nL_C)
        move_nmol(T_list,
                   src="complexU_N",\
                   dst="complexU_NPC",\
                   nmol=nU_N)
        move_nmol(T_list,
                   src="complexU_C",\
                   dst="complexU_NPC",\
                   nmol=nU_C)

    @register_update()
    def get_nmol_complex_NPC_to_complex_N_C(self, T_list):
        """
        Number of complexed cargo-importin released from the NPC to the nucleus and cytoplasm 
        (assumed 50-50 between nucleus and cytoplasm)

        Return: dictionary with number of molecules to add/subtract from each species
        """
        f= self.fraction_complex_NPC_to_complex_N_C_per_sec \
            * self.dt_sec # fractions are fine (conceptually, a random variable)
        nL= f * self.nmol["complexL_NPC"] 
        nU= f * self.nmol["complexU_NPC"] 
        move_nmol(T_list,
                   src="complexL_NPC",\
                   dst="complexL_N",\
                   nmol=0.5*nL)
        move_nmol(T_list,
                   src="complexL_NPC",\
                   dst="complexL_C",\
                   nmol=0.5*nL)
        move_nmol(T_list,
                   src="complexU_NPC",\
                   dst="complexU_N",\
                   nmol=0.5*nU)
        move_nmol(T_list,
                   src="complexU_NPC",\
                   dst="complexU_C",\
                   nmol=0.5*nU)

    @register_update()
    def get_nmol_cargo_bleached(self, T_list):
        '''
        Number of bleached molecules over time step in cytoplasm (both free and complexed)

        Return: dictionary with number of molecules to add/subtract from each species
        '''
        global N_A
        if self.sim_time_sec <= self.bleach_start_time_sec:
            return {}
        f= self.bleach_volume_L_per_sec \
            * N_A \
            * self.dt_sec
        c_freeL_C_M= self.get_concentration_M('freeL_C')
        n_free_C= f * c_freeL_C_M 
        c_complexL_C_M= self.get_concentration_M('complexL_C')
        n_complex_C= f * c_complexL_C_M 
        #print(f"Bleaching {n_free_C} free cargo molecules")
        #print(self.nmol["freeL_C"], f)
        assert(n_free_C <= self.nmol["freeL_C"])
        assert(n_complex_C <= self.nmol["complexL_C"])
        move_nmol(T_list,
                   src="freeL_C",\
                   dst="freeU_C",\
                   nmol=n_free_C)
        move_nmol(T_list,
                   src="complexL_C",\
                   dst="complexU_C",\
                   nmol=n_complex_C)

    ##########################
    # Individual update rules:
    ########################

    def get_nmol_T_summary(self, T_list):
        '''
        Summarize all transition by summing over a list of tuples of transitions

        @param T_list a list of tuples (src, dst, nmol), each representing a transfer of nmol 
                      molecules from a molecular species src to the molecular species dst
        @return a dictionary mapping from molecular species to total change in counts
        '''
        T= {}
        for src, dst, nmol in T_list:
            if src in T:
                T[src] -= nmol
            else:
                T[src] = -nmol

            if dst in T:
                T[dst] += nmol
            else:
                T[dst] = nmol

        return T

    def get_import_export_summary(self, T_list):
        '''
        Compute a summary of the gross number of cargo molecules imported and exported.
        Should be called BEFORE transitions are updated
        '''
        active_import = 0
        passive_import = 0
        active_export = 0
        passive_export = 0
        npc_docked_N = 0
        npc_docked_C = 0
        npc_undocked = 0

        for src, dst, nmol in T_list:
            # we aren't interested in Ran
            if "GTP" in src or "GDP" in src:
                continue

            # active import/export:
            if "NPC" in src:
                npc_undocked += nmol
                if "N" in dst:
                    active_import += (1-self.NPC_N2C_ratio) * nmol
                else:
                    active_export += self.NPC_N2C_ratio * nmol
            
            #passive export
            elif "N" in src and "C" in dst:
                passive_export += nmol
            #passive import
            elif "C" in src and "N" in dst:
                passive_import += nmol
            
            # for NPC_N2C_ratio
            if "NPC" in dst:
                if "N" in src:
                    npc_docked_N += nmol
                elif "C" in src:
                    npc_docked_C += nmol
        
        # update NPC_N2C_ratio
        total_NPC = sum([self.nmol[key] for key in self.nmol if "NPC" in key])
        total_NPC -= npc_undocked
        enumerator = self.NPC_N2C_ratio*total_NPC + npc_docked_N
        denominator = total_NPC + npc_docked_C + npc_docked_N
        self.NPC_N2C_ratio = enumerator/denominator
        
        total_N = sum([self.nmol[key] for key in self.nmol if "N" in key
                                                            and not "NPC" in key
                                                            and not "G" in key])
        total_C = sum([self.nmol[key] for key in self.nmol if "C" in key
                                                            and not "NPC" in key
                                                            and not "G" in key])
        self.nmol["import"] = (active_import + passive_import)/(total_N*self.dt_sec)
        self.nmol["export"] = (active_export + passive_export)/(total_C*self.dt_sec)


    def do_one_time_step(self):
        '''
        Update all state variables over a single time step
        '''
        # Compute transitions:
        T_list = []
        for update_rule in self._update_funcs:
            update_rule(self, T_list)
        T= self.get_nmol_T_summary(T_list)
        self.get_import_export_summary(T_list)

        # Update transitions:
        for key, value in T.items():
            if (key not in self.nmol):
                raise ValueError(f"can't update non-existent molecular species {key}")
            self.nmol[key] += value 
            if(self.nmol[key] < 0):
                if self.nmol[key]>-0.001:
                    T[key]= 0.0 
                else:
                    print(f"Negative key {key} value {self.nmol[key]} change {value}")
                    print(T)
                    print(self.nmol)
                    assert(self.nmol[key] >= 0)        
        # Update simulation clock:
        self.sim_time_sec += self.dt_sec

    def simulate(self, sim_time_sec, nskip_statistics= 10):
        ''' 
        simulate for approximately (and at least) sim_time_sec seconds
        @return actual time simulated
        '''
        # Computes number of steps and frames
        nsteps= int(np.ceil(sim_time_sec/self.dt_sec))
        nframes= ((nsteps-1)//nskip_statistics) + 1
        # Prepare statistics dictionary for all molecule types
        stats= { 'time_sec' : np.zeros(nframes) }
        for key in self.nmol.keys():
            stats[key]= np.zeros(nframes)
        for i in range(nsteps):
            self.do_one_time_step()
            if i % nskip_statistics == 0:
                si= i//nskip_statistics
                stats['time_sec'][si]= self.sim_time_sec
                for key, value in self.nmol.items():
                    stats[key][si]= value
        return stats

    ##########################
    # Debug utility functions
    ########################
    def get_total_RAN(self):
        RAN = self.nmol["GDP_C"] + self.nmol["GTP_C"] + self.nmol["GDP_N"] + self.nmol["GTP_N"]
        return RAN
    def get_total_cargoL_nmol(self):
        return self.nmol["complexL_C"] + self.nmol["freeL_C"] + \
               self.nmol["complexL_N"] + self.nmol["freeL_N"] + \
               self.nmol["complexL_NPC"] 
    def get_total_cargoU_nmol(self):
        return self.nmol["complexU_C"] + self.nmol["freeU_C"] + \
               self.nmol["complexU_N"] + self.nmol["freeU_N"] + \
               self.nmol["complexU_NPC"]
    def get_total_cargo_nmol(self):
        return self.get_total_cargoL_nmol() + self.get_total_cargoU_nmol()

