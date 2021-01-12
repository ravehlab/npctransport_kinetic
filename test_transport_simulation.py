# Synopsis: unit tests for TransportSimulation

import unittest
import transport_simulation 
import numpy as np
from numpy import array

# +
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
            {'time_sec': array([1.000e-03, 2.501e+00, 5.001e+00, 7.501e+00]), 
             'complexL_NPC_N_import': array([    0.        , 15239.18509368, 19158.28537664, 20402.53909047]), 
             'complexL_NPC_C_import': array([    0.        , 15287.69916179, 19211.99563442, 20456.99930688]), 
             'complexU_NPC_N_import': array([0., 0., 0., 0.]), 
             'complexU_NPC_C_import': array([0., 0., 0., 0.]), 
             'complexL_NPC_N_export': array([    0.        ,  6490.57957647, 12676.39355691, 16353.10477062]), 
             'complexL_NPC_C_export': array([    0.        ,  6471.41018374, 12647.33483225, 16318.21182419]), 
             'complexU_NPC_N_export': array([0., 0., 0., 0.]), 
             'complexU_NPC_C_export': array([0., 0., 0., 0.]), 
             'complexL_C': array([   30.11      ,  9006.61140819, 15560.61813652, 21813.31804998]), 
             'freeL_C': array([301069.889     , 235126.55119532, 184485.65932088, 145763.83104942]), 
             'complexU_C': array([0., 0., 0., 0.]), 
             'freeU_C': array([0., 0., 0., 0.]), 
             'complexL_N': array([   0.        , 1398.96732942, 3346.73873462, 5487.54771013]), 
             'freeL_N': array([1.00000000e-03, 1.20789961e+04, 3.40129744e+04, 5.45044482e+04]), 
             'complexU_N': array([0., 0., 0., 0.]), 
             'freeU_N': array([0., 0., 0., 0.]), 
             'import_L': array([0.        , 1.59948125, 0.69655085, 0.44976187]), 
             'export_L': array([0.        , 0.02649327, 0.06320976, 0.09736769]), 
             'import_U': array([0., 0., 0., 0.]), 
             'export_U': array([0., 0., 0., 0.]), 
             'GTP_N': array([78196.06103896, 74770.68735182, 64671.25396217, 58224.01573541]), 
             'GTP_C': array([50.83506494, 40.991402  , 43.39131566, 42.92777012]), 
             'GDP_N': array([ 93.83370909, 156.21858551, 156.21617624, 156.2166379 ]), 
             'GDP_C': array([78231.27018701, 81604.10266066, 91701.13854594, 98148.83985657])
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
            

# -



if __name__ == '__main__':
    unittest.main()
\
