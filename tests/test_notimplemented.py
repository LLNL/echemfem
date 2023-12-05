import pytest
from firedrake import *
from echemfem import EchemSolver

flows = [["electroneutrality"],
         ["poisson", "migration", "electroneutrality"],
         ["poisson", "migration", "electroneutrality full"],
         ["migration"]
         ]


class Solver(EchemSolver):
    def __init__(self, flow, family="CG"):
        mesh = UnitIntervalMesh(2)
        conc_params = []
        physical_params = []
        physical_params = {"flow": flow,
                           }

        super().__init__(conc_params, physical_params, mesh, family=family)

    def set_boundary_markers(self):
        self.boundary_markers = {}


def _test_flow(flow):
    # with pytest.raises(NotImplementedError):
    #    solver = Solver(flow)
    solver = Solver(flow)


def test_flows():
    for flow in flows:
        try:
            _test_flow(flow)
        except NotImplementedError as e:
            print(e)


test_flows()
