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
        ts.bleach_start_time_sec= 100.0
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

if __name__ == '__main__':
    unittest.main()
