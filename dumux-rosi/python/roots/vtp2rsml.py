import sys; sys.path.append("../modules/"); sys.path.append("../../../CPlantBox/");  sys.path.append("../../../CPlantBox/src/python_modules")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml_writer as rsml
import vtk_tools as vt
import vtk_plot as vp

from math import *
import numpy as np
import matplotlib.pyplot as plt

""" 
Converts a DuMux output vtp to a RSML
"""

file_in = "test.vtp"  # ../../grids/RootSystem8.vtp"
file_out = "test.rsml"  # "../../grids/RootSystem8.rsml"

""" read vtp """
pd = vt.read_vtp(file_in)

meta = rsml.Metadata()
meta.unit = "m"
meta.add_property(rsml.Property("radius [m]", "float", "m", None))

order_id = 4

vt.write_rsml(file_out, pd, order_id, meta)  # meta is optional now

print("fin")
