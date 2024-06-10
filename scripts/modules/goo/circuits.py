class Gene:
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return self.name


class Circuit:
    def __init__(self, *components):
        self.circuit = "\n".join(components)

    def __str__(self):
        return self.circuit


def deg_first_order(x: Gene, k: float):
    """d[x]/dt = -k[x]"""
    return f"{x} -> ; {k} * {x}"


def prod_and(
    z: Gene,
    x: Gene,
    y: Gene,
    k: float,
    nx: float = 2,
    ny: float = 2,
    s: Gene = "",
    a0: float = 0,
):
    """d[z]/dt = k([x] AND [y]), optional substrate s consumed, and optional leaky factor a0."""
    return f"{s} -> {z}; {a0} + {k} * {x}^{nx} * {y}^{ny} / (1 + {x}^{nx}) / (1 + {y}^{ny})"


def prod_activation(
    y: Gene,
    x: Gene,
    kcat: float,
    Km: float = 1,
    n: float = 2,
    s: Gene = "",
    a0: float = 0,
):
    """d[y]/dt = kcat * [x]**n / (Km + [x]**n), optional substrate s consumed, and optional leaky factor a0."""
    return f"{s} -> {y}; {a0} + {kcat} * {x}^{n} / ({Km} + {x}^{n})"


def prod_repression(
    y: Gene,
    x: Gene,
    kcat: float,
    Km: float = 1,
    n: float = 2,
    s: Gene = "",
    a0: float = 0,
):
    """d[y]/dt = kcat / (Km + [x]**n), optional substrate s consumed, and optional leaky factor a0."""
    return f"{s} -> {y}; {a0} + {kcat} / ({Km} + {x}^{n})"
