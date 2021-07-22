# Synopsis: unit tests for TransportSimulation
import numpy as np
import os
import pickle
import tempfile
import stats_grid
import subprocess
import unittest
import transport_simulation 

class TestStatsGrid(unittest.TestCase):
     
    def test_simple_grid(self):
        print("Testing simple grid")
        pkl_file = tempfile.NamedTemporaryFile(suffix=".pkl")
        print(f"Pickle file {pkl_file.name}")
        stats_grid.get_stats_on_grid(
            output="test",
            passive_range=(0.01,0.15),
            npc_traverse_range=(1,1000),
            k_on_range=(0.001,5),
            nx=2,
            ny=2,
            n_passive=2,
            cargo_concentration_M=0.1e-6,
            Ran_concentration_M=20e-6,
            v_N_L=627e-15,
            v_C_L=2194e-15,
            equilibration_time_sec = 0.2,            
            pickle_file = pkl_file.name,
            number_of_processors=4)
        with open(pkl_file.name,"rb") as F:
            output = pickle.load(F)
            for i, oi in enumerate(output):
                print(f"i={i}")
                print(oi)
            stats_by_passive, TSs, tags = output 
            self.assertAlmostEqual(stats_by_passive[0.01]['complexL_C'][0][0], 1.23323586e4, places=0)
            self.assertAlmostEqual(stats_by_passive[0.01]['complexL_C'][0][1], 8.97198131e3, places=0)
            self.assertAlmostEqual(stats_by_passive[0.01]['complexL_C'][1][0], 4.00061993e+07, places=0)
            self.assertAlmostEqual(stats_by_passive[0.01]['complexL_C'][1][1],  3.06235254e+07, places=0)
            self.assertEqual(len(TSs), 2)
            self.assertAlmostEqual(TSs[0.01].v_N_L, 627e-15)
            self.assertAlmostEqual(TSs[0.01].v_C_L, 2194e-15)



if __name__ == '__main__':
    unittest.main()


