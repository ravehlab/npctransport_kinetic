# Synopsis: unit tests for TransportSimulation

import unittest
import transport_simulation 
import numpy as np
from numpy import array

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
            {'time_sec': array([1.000e-03, 2.501e+00, 5.001e+00, 7.501e+00]), 'complexL_NPC_N_import': array([    0.        , 26716.61622277, 27661.69858981, 26739.33829553]), 'complexL_NPC_C_import': array([    0.        , 26841.59293957, 27778.30377937, 26846.93361103]), 'complexU_NPC_N_import': array([0., 0., 0., 0.]), 'complexU_NPC_C_import': array([0., 0., 0., 0.]), 'complexL_NPC_N_export': array([    0.        ,  8730.19487584, 14009.35264709, 17086.61559209]), 'complexL_NPC_C_export': array([    0.        ,  8706.54011213, 13978.39393505, 17050.47109701]), 'complexU_NPC_N_export': array([0., 0., 0., 0.]), 'complexU_NPC_C_export': array([0., 0., 0., 0.]), 'complexL_C': array([1.68164350e+02, 2.10961044e+05, 3.48534180e+05, 4.39349880e+05]), 'freeL_C': array([1681475.33465   , 1323004.67669619, 1062099.67098618,
        871297.61991343]), 'complexU_C': array([0., 0., 0., 0.]), 'freeU_C': array([0., 0., 0., 0.]), 'complexL_N': array([    0.        ,  6292.36410467, 14737.13433589, 22944.203708  ]), 'freeL_N': array([1.00000000e-03, 7.03904708e+04, 1.72844766e+05, 2.60328437e+05]), 'complexU_N': array([0., 0., 0., 0.]), 'freeU_N': array([0., 0., 0., 0.]), 'nuclear_importL_per_sec': array([7.63484083e-06, 5.09284127e-01, 5.32553385e-01, 5.31471750e-01]), 'nuclear_exportL_per_sec': array([0.        , 0.11357389, 0.07453442, 0.06020191]), 'nuclear_importU_per_sec': array([0., 0., 0., 0.]), 'nuclear_exportU_per_sec': array([0., 0., 0., 0.]), 'GTP_N': array([362107.91342657, 334681.57208627, 292146.25547386, 269547.11807043]), 'GTP_C': array([235.40545455, 199.01188995, 197.02947705, 194.66540413]), 'GDP_N': array([434.52225287, 723.40299972, 723.40495259, 723.40731827]), 'GDP_C': array([362270.95886601, 389444.81302406, 431982.1100965 , 454583.60920717])}
        ts= transport_simulation.TransportSimulation()
        ts.set_params(fraction_complex_NPC_to_free_N_per_M_GTP_per_sec= 0.010e+6)
        ts.bleach_start_time_sec= 100.0
        ts.dt_sec= 1e-3 
        sim_flags= dict()#rate_free_to_complex_per_sec=1.0,
           #max_passive_diffusion_rate_nmol_per_sec_per_M=2e7)
        sim_time_sec= 10.0
        stats=ts.simulate(sim_time_sec, nskip_statistics= 2500)
#        print(stats)
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
        
    def test_get_compartment(self):
        ts= transport_simulation.TransportSimulation()
        self.assertEqual(ts.get_compartment('complexL_NPC_N_import'), "NPC")
        self.assertEqual(ts.get_compartment('complexL_NPC_C_import'), "NPC")
        self.assertEqual(ts.get_compartment('complexU_NPC_N_export'), "NPC")
        self.assertEqual(ts.get_compartment('complexL_NPC_C_export'), "NPC")
        self.assertEqual(ts.get_compartment('complexL_C'), "C")
        self.assertEqual(ts.get_compartment('freeL_C'), "C")
        self.assertEqual(ts.get_compartment('complexU_N'), "N")
        self.assertEqual(ts.get_compartment('freeL_N'), "N")
        self.assertEqual(ts.get_compartment('GTP_N'), "N")
        self.assertEqual(ts.get_compartment('GDP_N'), "N")
        self.assertEqual(ts.get_compartment('GTP_C'), "C")
        self.assertEqual(ts.get_compartment('GDP_C'), "C")
        
    def test_set_N_C_volume(self):
        ts= transport_simulation.TransportSimulation()
        n1= ts.get_nmol('freeL_C')
        C1= ts.get_concentration_M('freeL_C')
        print(ts.get_v_C_L())
        ts.set_v_C_L(1000, fix_concentration= True)
        print(ts.get_v_C_L())
        n2= ts.get_nmol('freeL_C')
        C2= ts.get_concentration_M('freeL_C')
        print(f"n1={n1}, n2={n2}, C1={C1}, C2={C2}")
        self.assertAlmostEqual(C1, C2)



if __name__ == '__main__':
    unittest.main()


