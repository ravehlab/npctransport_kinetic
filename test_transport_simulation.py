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
            {'time_sec': [1.000e-03, 2.501e+00, 5.001e+00, 7.501e+00],
             'complexL_NPC_N_import': [    0.        , 21623.81698188, 21278.51029702, 21176.81632042],
             'complexL_NPC_C_import': [  625.        , 21681.06379079, 21334.39189677, 21231.62598839],
             'complexU_NPC_N_import': [0., 0., 0., 0.],
             'complexU_NPC_C_import': [0., 0., 0., 0.],
             'complexL_NPC_N_export': [    0.        , 13545.21145814, 17168.79328743, 19146.95996426],
             'complexL_NPC_C_export': [    0.        , 13513.62995572, 17132.5458864 , 19107.47897038],
             'complexU_NPC_N_export': [0., 0., 0., 0.],
             'complexU_NPC_C_export': [0., 0., 0., 0.],
             'complexL_C': [74668.81875   , 20207.78255697, 25699.87742521, 30248.80231033],
             'freeL_C': [225806.1805    , 178129.82134661, 141282.393471  , 113145.48187596],
             'complexU_C': [0., 0., 0., 0.],
             'freeU_C': [0., 0., 0., 0.],
             'complexL_N': [   0.        , 4433.24433551, 6522.12774267, 8448.79229082],
             'freeL_N': [7.50000000e-04, 2.79654296e+04, 5.06813600e+04, 6.85940423e+04],
             'complexU_N': [0., 0., 0., 0.],
             'freeU_N': [0., 0., 0., 0.],
             'import_L': [0.        , 0.90171421, 0.48877104, 0.35666296],
             'export_L': [0.        , 0.06811956, 0.10259178, 0.13324456],
             'import_U': [0., 0., 0., 0.],
             'export_U': [0., 0., 0., 0.],
             'GTP_N': [78196.06103896, 63317.01960248, 56680.03947056, 53731.94361351],
             'GTP_C': [50.83506494, 45.2034439 , 43.17588964, 42.61883649],
             'GDP_N': [ 93.83370909, 156.21436088, 156.21638929, 156.21694678],
             'GDP_C': [ 78231.27018701,  93053.56259274,  99692.56825052, 102641.22060322]
             }
        ts= transport_simulation.TransportSimulation()
        ts.set_params(fraction_complex_NPC_to_free_N_per_M_GTP_per_sec= 0.010e+6)
        ts.bleach_start_time_sec= 100.0
        ts.dt_sec= 1e-3 
        sim_flags= dict()#rate_free_to_complex_per_sec=1.0,
           #max_passive_diffusion_rate_nmol_per_sec_per_M=2e7)
        sim_time_sec= 10.0
        stats=ts.simulate(sim_time_sec, nskip_statistics= 2500)
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
