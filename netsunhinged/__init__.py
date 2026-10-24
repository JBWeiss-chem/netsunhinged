"""NetsUnhinged - large-amplitude flexibility of framework materials (MOFs, ZIFs, zeolites)
from rigid-body / hinge models.

Implements and generalises the approach of Sarkisov, Martin, Haranczyk & Smit,
"On the Flexibility of Metal-Organic Frameworks", J. Am. Chem. Soc. 2014, 136, 2228.
"""

from .cif import read_cif
from .structure import Structure
from .model import MechanicalModel, ModelOptions, preset

__all__ = ["read_cif", "Structure", "MechanicalModel", "ModelOptions", "preset", "load"]
__version__ = "0.1.0"


def load(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        return Structure(read_cif(fh.read()))
