# Synopsis: unit tests for TransportSimulation

import unittest
import transport_simulation 
import numpy as np

class TestTransport_Simulation(unittest.TestCase):
     
    def test_fL_to_L(self):
        self.assertAlmostEqual(2e-15,
                               transport_simulation.fL_to_L(2.0),
                               delta=1e-18)

    def test_simple_simulation(self):
        # TODO: make this more interesting
        ts= transport_simulation.TransportSimulation()
        ts.dt_sec= 1e-3
        sim_time_sec= 0.1
        n_skip= 1
        n= int(np.ceil(sim_time_sec / ts.dt_sec))
        attributes = ['GDP_N', 'GDP_C', 'GTP_N', 'GTP_C']
        print(f"Running for n={n} time steps")
        n_frames= int(np.floor(n/n_skip)) + 1 
        c_arr = np.zeros((len(attributes), n_frames))
        i_frame= 0
        print("Cargo before: {}".format(ts.get_total_cargo_nmol()))
        for i in range(n):
            ts.do_one_time_step()
            if (i % n_skip) == 0:
                for j, attr in enumerate(attributes):
                    c_arr[j,i_frame] = ts.get_nmol(attr)
                i_frame += 1
        print("Cargo after: {}".format(ts.get_total_cargo_nmol()))
        #return c_arr[:,:-1]

    def test_simulation(self):
        
        expected_stats= \
            {'time_sec': np.array([1.000e-03, 2.501e+00, 5.001e+00, 7.501e+00]), 
             'complexL_NPC_N_import': np.array([    0.        , 21624.82871907, 21267.58400813, 21155.07290699]), 
             'complexL_NPC_C_import': np.array([  625.        , 21682.07622186, 21323.44677205, 21209.84595497]), 
             'complexU_NPC_N_import': np.array([0., 0., 0., 0.]), 
             'complexU_NPC_C_import': np.array([0., 0., 0., 0.]), 
             'complexL_NPC_N_export': np.array([    0.        , 13466.96758474, 17038.67805727, 18980.21813115]), 
             'complexL_NPC_C_export': np.array([    0.        , 13435.5924849 , 17002.72500676, 18941.09400388]), 
             'complexU_NPC_N_export': np.array([0., 0., 0., 0.]), 
             'complexU_NPC_C_export': np.array([0., 0., 0., 0.]), 
             'complexL_C': np.array([74668.81875   , 20101.80884471, 25377.20385103, 29643.56416695]), 
             'freeL_C': np.array([225806.1805    , 178125.28595829, 141255.45742185, 113072.97871725]), 
             'complexU_C': np.array([0., 0., 0., 0.]), 
             'freeU_C': np.array([0., 0., 0., 0.]), 
             'complexL_N': np.array([   0.        , 4382.17353   , 6392.74622079, 8215.18778046]), 
             'freeL_N': np.array([7.50000000e-04, 2.82812667e+04, 5.14421587e+04, 6.98820383e+04]), 
             'complexU_N': np.array([0., 0., 0., 0.]), 
             'freeU_N': np.array([0., 0., 0., 0.]), 
             'import_L': np.array([0.        , 0.89455947, 0.48334772, 0.35167244]), 
             'export_L': np.array([0.        , 0.06776401, 0.10202806, 0.13271156]), 
             'import_U': np.array([0., 0., 0., 0.]), 
             'export_U': np.array([0., 0., 0., 0.]), 
             'GTP_N': np.array([78196.06103896, 63347.38643075, 56758.89447587, 53855.74247128]), 
             'GTP_C': np.array([50.83506494, 45.15307355, 43.10669637, 42.53502578]), 
             'GDP_N': np.array([ 93.83370909, 156.2144111 , 156.21645831, 156.2170304 ]), 
             'GDP_C': np.array([ 78231.27018701,  93023.2460846 ,  99613.78236945, 102517.50547255])
            }
        ts= transport_simulation.TransportSimulation()
        ts.set_params(fraction_complex_NPC_to_free_N_per_M_GTP_per_sec= 0.010e+6)
        ts.bleach_start_time_sec= 100.0
        ts.dt_sec= 1e-3 
        sim_flags= dict()#rate_free_to_complex_per_sec=1.0,
           #max_passive_diffusion_rate_nmol_per_sec_per_M=2e7)
        sim_time_sec= 10.0
        stats=ts.simulate(sim_time_sec, nskip_statistics= 2500)
        print(stats)
        for key,expected_v in expected_stats.items():
            for x,y in zip(stats[key], expected_v):
                self.assertTrue(np.isclose(x, y, atol= 1e-12, rtol=1e-5))

    def test_constant_cargo(self, nsteps=1000):
        ts= transport_simulation.TransportSimulation()
        for i in range(nsteps):
            RAN = ts.get_total_RAN()
            cargo = ts.get_total_cargo_nmol()
            ts.do_one_time_step()
            RAN_after = ts.get_total_RAN()
            cargo_after = ts.get_total_cargo_nmol()
            self.assertAlmostEqual(RAN, RAN_after, places=1)
            self.assertAlmostEqual(cargo, cargo_after, places=1)



if __name__ == '__main__':
    unittest.main()
\
