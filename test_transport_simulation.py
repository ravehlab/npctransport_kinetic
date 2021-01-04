# Synopsis: unit tests for TransportSimulation

import unittest
import transport_simulation as ts

class TestTransport_Simulation(unittest.TestCase):

    def setUp(self):
        self.sim= ts.TransportSimulation()

    def _get_default_simulation():
        return self.sim
        
    def test_fL_to_L(self):
        self.assertAlmostEqual(2e-15,
                               ts.fL_to_L(2.0),
                               delta=1e-18)
        

if __name__ == '__main__':
    unittest.main()
