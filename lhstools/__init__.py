# from .benchmark import Benchmark
# from .fapar import FAPAR
# from .emdi import EMDI
# from .uptake import Uptake
# from .globalveg import GlobalVeg
# from .globalsoil import GlobalSoil
# from .fluxnet import FLUXNET
# from .uptake80 import Uptake80
# from .uptake90 import Uptake90
# from .uptake00 import Uptake00
# from .uptake02 import Uptake02
# from .uptaketot import UptakeTot
# from .soilcmap import SoilCMap
# from .soilcmaplow import SoilCMapLow
# from .totcmap import TotCMap
# from .emdi import EMDI

# Allow from lhstools import *
# from os.path import dirname, basename, isfile
# import glob
# modules = glob.glob(dirname(__file__)+"/*.py")
# __all__ = [ basename(f)[:-3] for f in modules if isfile(f)]
import os, sys

path = os.path.dirname(os.path.abspath(__file__))

for py in [f[:-3] for f in os.listdir(path) if f.endswith('.py') and f != '__init__.py']:
    mod = __import__('.'.join([__name__, py]), fromlist=[py])
    classes = [getattr(mod, x) for x in dir(mod) if isinstance(getattr(mod, x), type)]
    for cls in classes:
        setattr(sys.modules[__name__], cls.__name__, cls)
