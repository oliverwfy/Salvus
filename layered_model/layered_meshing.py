import numpy as np
import xarray as xr

import salvus.namespace as sn
import salvus.mesh.layered_meshing as lm

from Salvus.material_ela_constants.Elastic_Material import Austenite



material = Austenite()

model = lm.material.elastic.hexagonal.TensorComponents(*material.params)
d_3d = sn.domain.dim3.BoxDomain(x0=0, x1=1, y0=0, y1=1, z0=0, z1=1)

