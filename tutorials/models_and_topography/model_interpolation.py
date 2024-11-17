import os
from functools import partial
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

import salvus.namespace as sn
import salvus.toolbox.toolbox as st
from my_mesh import *

# Salvus site name
SALVUS_FLOW_SITE_NAME = 'oliver_wsl'
ranks = 4


# Directories in WSL
PROJECT_DIR = '/home/oliver/workspace/Salvus/tutorials/Project/model_interpolation'

# material properties (linear isotropic elastic solid)
matls = Isotropic_Material(RHO=2200.0, VP=6000.0, VS=3000.0)

# 2D box domain parameters (length in m)
x_length = 10 * 1e-3
y_length = 10 * 1e-3

# centre frequency 
f_c = 5*1e6


tensor_order = 4                        # Order of the GLL model and shape mapping.
elements_per_wavelength = 2             # Elements per mimimum wavelength.
abs_bdry = None

m = get_basic_mesh(matlt = matls, x_max=x_length, y_max=y_length
                   , centre_freq=f_c, elements_per_wavelength = elements_per_wavelength, 
                   tensor_order=tensor_order, abs_bdry=None)




# shape (m x n x d)  m=#ele, n=(tensor_order+1)^dim d = dim
mess_nodes = m.get_element_nodes()       
m.get_element_nodes().shape



# replace existing elemental_fields with new material properties
m.elemental_fields.clear()
rho=2200.0 
vp=6000.0 
vs=3000.0


# define Lamé constants by wave velocities 
mu = rho * vs**2
lam = rho * vp**2 - 2 * mu

# Attach parameter to the nodes of each element.
par_template = np.ones_like(m.get_element_nodes()[:, :, 0])
m.attach_field("LAMBDA", par_template * lam)
m.attach_field("MU", par_template * mu)
m.attach_field("RHO", par_template * rho)
m.attach_field("fluid", np.zeros(m.nelem))



# set input time function
stf = sn.simple_config.stf.Ricker(center_frequency=f_c)
src_vector_2d = sn.simple_config.source.cartesian.VectorPoint2D(
    fx=0, fy=1, x=0, y=0, source_time_function=stf
)


# set up waveform 
w = sn.simple_config.simulation.Waveform()
w.domain.dimension = 2
w.output.volume_data.format = "hdf5"
w.output.volume_data.filename = "output.h5"
w.output.volume_data.sampling_interval_in_time_steps = 100


# Attach the mesh and set some custom output.
w.set_mesh(m)
w.output.volume_data.fields = ["displacement"]


w.physics.wave_equation.point_source = [src_vector_2d]

# Run the solver.
output_folder = Path("./parameterization/elastic_lambdamurho")
output_file = output_folder / "output.h5"


run_salvus = partial(
    sn.api.run, ranks=ranks, get_all=True, site_name=SALVUS_FLOW_SITE_NAME
)

run_salvus(input_file=w, output_folder=output_folder)

# Visualize the results.
f, ax = plt.subplots(1, 1)
ax.set_aspect("equal")
t, de1 = st.visualize_wavefield_2d(output_file, "displacement")
ax.tricontourf(t, de1[-1, :])