from goo.circuits import *
from math import e
import unittest

x = Gene("x")
y = Gene("y")
z = Gene("z")
s = Gene("s")
engine = TelluriumEngine()


class TestCircuits(unittest.TestCase):
    def concs_almost_equal(self, expected, actual, places=4):
        for gene, conc in expected.items():
            self.assertAlmostEqual(conc, actual[str(gene)], places=places)

    def test_dfo(self):
        engine.load_circuits(DegFirstOrder(x, 0.2))
        concs = {x: 2}

        engine.step(concs)
        self.concs_almost_equal({x: 1.63746}, engine.retrieve_concs())

    def test_activation(self):
        engine = TelluriumEngine()
        engine.load_circuits(ProdActivation(y, x, kcat=0.2, Km=2, s=s))
        concs = {y: 2, x: 5, s: 10}

        engine.step(concs)
        self.concs_almost_equal({y: 2.185185, s: 9.8148}, engine.retrieve_concs())

    def test_dfo_activation(self):
        engine.load_circuits(
            DegFirstOrder(x, 0.2), ProdActivation(y, x, kcat=0.2, Km=2, s=s)
        )
        concs = {y: 2, x: 5, s: 10}
        engine.step(concs)
        self.concs_almost_equal(
            {x: 4.09365, y: 2.1821, s: 9.81789}, engine.retrieve_concs()
        )


if __name__ == "__main__":
    unittest.main()
