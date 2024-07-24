import antimony
import roadrunner
from dataclasses import dataclass

from .molecule import DiffusionSystem


class Gene:
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return self.name


class Circuit:
    pass


@dataclass
class DegFirstOrder(Circuit):
    x: Gene
    k: float
    """d[x]/dt = -k[x]"""


@dataclass
class ProdAnd(Circuit):
    """d[z]/dt = k([x] AND [y]), optional substrate s consumed, and optional leaky factor a0."""

    z: Gene
    x: Gene
    y: Gene
    k: float
    nx: float = 2
    ny: float = 2
    s: Gene = None
    a0: float = 0


@dataclass
class ProdActivation(Circuit):
    """d[y]/dt = kcat * [x]**n / (Km + [x]**n), optional substrate s consumed, and optional leaky factor a0."""

    y: Gene
    x: Gene
    kcat: float
    Km: float = 1
    n: float = 2
    s: Gene = None
    a0: float = 0


@dataclass
class ProdRepression(Circuit):
    """d[y]/dt = kcat / (Km + [x]**n), optional substrate s consumed, and optional leaky factor a0."""

    y: Gene
    x: Gene
    kcat: float
    Km: float = 1
    n: float = 2
    s: Gene = None
    a0: float = 0


class CircuitEngine:
    def __init__(self):
        pass

    def load_circuits(self, *circuits: Circuit):
        """Load engine with model network."""
        pass

    def step(self, metabolite_concs: dict):
        """Take 1 frame step of circuit calculations to calculate new metabolite concentrations."""
        pass

    def retrieve_concs(self):
        """Get current metabolite concentrations."""
        pass


class TelluriumEngine(CircuitEngine):
    def __init__(self):
        self.model = None
        self.result = None

    def load_circuits(self, *circuits):
        model = []
        for circuit in circuits:
            model.append(self._encode(circuit))
        self.model = "\n".join(model)

    @staticmethod
    def _encode(circuit):
        match circuit:
            case DegFirstOrder(x, k):
                return f"{x} -> ; {k} * {x}"
            case ProdAnd(z, x, y, k, nx, ny, s, a0):
                s = "" if s is None else s
                return f"{s} -> {z}; {a0} + {k} * {x}^{nx} * {y}^{ny} / (1 + {x}^{nx}) / (1 + {y}^{ny})"
            case ProdActivation(y, x, kcat, Km, n, s, a0):
                s = "" if s is None else s
                return f"{s} -> {y}; {a0} + {kcat} * {x}^{n} / ({Km} + {x}^{n})"
            case ProdRepression(y, x, kcat, Km, n, s, a0):
                s = "" if s is None else s
                return f"{s} -> {y}; {a0} + {kcat} / ({Km} + {x}^{n})"

    def load_sbml(self, sbml):
        pass

    def step(self, metabolite_concs, iter=5):
        if self.model is None:
            return
        model_full = self._get_model_full(metabolite_concs)

        # Replacement of Tellurium with lower-level functions
        antimony.clearPreviousLoads()
        antimony.freeAll()
        code = antimony.loadAntimonyString(model_full)
        mid = antimony.getMainModuleName()
        rr = roadrunner.RoadRunner(antimony.getSBMLString(mid))

        result = rr.simulate(0, 1, iter)
        self.result = result

    def _get_model_full(self, metabolite_concs):
        prefix = "model cell"
        suffix = "end"
        gene_levels = "\n".join(
            [
                f"{metabolite} = {level}"
                for metabolite, level in metabolite_concs.items()
            ]
        )
        return "\n".join([prefix, self.model, gene_levels, suffix])

    def retrieve_concs(self):
        if self.result is None:
            return {}
        colnames = map(lambda name: name[1:-1], self.result.colnames[1:])
        new_concs = dict(zip(colnames, self.result[-1][1:]))
        return new_concs


class GeneRegulatoryNetwork:
    def __init__(self):
        self.metabolite_concs = {}
        self.diffusion_system: DiffusionSystem = None
        self.circuit_engine: CircuitEngine = TelluriumEngine()

    def copy(self):
        pass

    @property
    def concs(self):
        return self.metabolite_concs

    @concs.setter
    def concs(self, concs):
        self.metabolite_concs = concs

    def set_diffusion_system(self, system):
        self.diffusion_system = system

    def load_circuits(self, *circuits):
        self.circuit_engine.load_circuits(*circuits)

    def update_concs(self, center=None, radius=None):
        """Update network concentrations based on underlying diffusion systems, cell center and radius."""
        if self.diffusion_system is not None:
            if center is None or radius is None:
                raise ValueError(
                    "If there is a diffusion system, center or radius must be supplied as arguments."
                )
            self.update_signaling_concs(center, radius)
        self.update_metabolite_concs()

    def update_metabolite_concs(self, iter=5):
        self.circuit_engine.step(self.metabolite_concs, iter=iter)
        new_concs = self.circuit_engine.retrieve_concs()
        self.metabolite_concs.update(new_concs)

    def update_signaling_concs(self, center, radius):
        signaling_concs = self.diffusion_system.get_all_concentrations(center, radius)
        self.metabolite_concs.update(signaling_concs)
