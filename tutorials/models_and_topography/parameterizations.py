import os
from functools import partial
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

import salvus.namespace as sn
import salvus.toolbox.toolbox as st

SALVUS_FLOW_SITE_NAME = 'oliver_mac'
PROJECT_DIR = "/Users/Oliver/Library/CloudStorage/OneDrive-UniversityofBristol/Salvus_project/model_and_topography"




def get_basic_mesh(dim: int, epd: int = 20) -> sn.UnstructuredMesh:
    """Get a simple mesh to outline allowed parameter types.

    Parameters
    ----------
    dim : int
        Dimension of the mesh.
    epd : int, optional
        Elements per dimension, by default 20

    Returns
    -------
    um.UnstructuredMesh
        An unstructured mesh free of any parameters.

    """
    x = 2.0
    h = x / float(epd)
    vp = 1000.0

    if dim == 2:
        mesh = sn.simple_mesh.basic_mesh.CartesianHomogeneousAcoustic2D(
            x_max=x,
            y_max=x,
            rho=1000.0,
            vp=vp,
            max_frequency=0.5 * vp / h,
            tensor_order=4,
        ).create_mesh()

    elif dim == 3:
        mesh = sn.simple_mesh.basic_mesh.CartesianHomogeneousAcoustic3D(
            x_max=x,
            y_max=x,
            z_max=x,
            rho=1000.0,
            vp=vp,
            max_frequency=0.5 * vp / h,
            tensor_order=4,
        ).create_mesh()

    # Delete the material properties as they will be added later on.
    mesh.elemental_fields.clear()
    return mesh

vs = 500.0
vp = 1000.0
rho = 1000.0


# define a set of sources
stf = sn.simple_config.stf.Ricker(center_frequency=2e3)

src_scalar_2d = sn.simple_config.source.cartesian.ScalarPoint2D(
    f=1, x=1, y=1, source_time_function=stf
)
src_scalar_3d = sn.simple_config.source.cartesian.ScalarPoint3D(
    f=1, x=1, y=1, z=1, source_time_function=stf
)
src_vector_2d = sn.simple_config.source.cartesian.VectorPoint2D(
    fx=1, fy=1, x=1, y=1, source_time_function=stf
)
src_vector_3d = sn.simple_config.source.cartesian.VectorPoint3D(
    fx=1, fy=1, fz=1, x=1, y=1, z=1, source_time_function=stf
)

run_salvus = partial(
    sn.api.run, ranks=4, get_all=True, site_name=SALVUS_FLOW_SITE_NAME
)

w = sn.simple_config.simulation.Waveform()

w.domain.dimension = 2
w.output.volume_data.format = "hdf5"
w.output.volume_data.filename = "output.h5"
w.output.volume_data.sampling_interval_in_time_steps = 100


# Generate the mesh.
# Velocities and density
m = get_basic_mesh(2)

# Attach parameter to the nodes of each element.
par_template = np.ones_like(m.get_element_nodes()[:, :, 0])
m.attach_field("VP", par_template * vp)
m.attach_field("VS", par_template * vs)
m.attach_field("RHO", par_template * rho)
m.attach_field("fluid", np.zeros(m.nelem))

# Attach the mesh and set some custom output.
w.set_mesh(m)
w.output.volume_data.fields = ["displacement"]
w.physics.wave_equation.point_source = [src_vector_2d]

# Run the solver.
output_folder = Path(PROJECT_DIR , "elastic_vpvsrho")
output_file = output_folder / "output.h5"
run_salvus(input_file=w, output_folder=output_folder)

# Visualize the results.
f, ax = plt.subplots(1, 1)
ax.set_aspect("equal")
t, de0 = st.visualize_wavefield_2d(output_file, "displacement")
ax.tricontourf(t, de0[-1, :])

# Generate the mesh.
# Lame parameters and density
m = get_basic_mesh(2)

# Attach parameter to the nodes of each element.
mu = rho * vs**2
lam = rho * vp**2 - 2 * mu
par_template = np.ones_like(m.get_element_nodes()[:, :, 0])
m.attach_field("LAMBDA", par_template * lam)
m.attach_field("MU", par_template * mu)
m.attach_field("RHO", par_template * rho)
m.attach_field("fluid", np.zeros(m.nelem))

# Attach the mesh and set some custom output.
w.set_mesh(m)
w.output.volume_data.fields = ["displacement"]
w.physics.wave_equation.point_source = [src_vector_2d]

# Run the solver.
output_folder = Path(PROJECT_DIR , "elastic_lambdamurho")
output_file = output_folder / "output.h5"
run_salvus(input_file=w, output_folder=output_folder)

# Visualize the results.
f, ax = plt.subplots(1, 1)
ax.set_aspect("equal")
t, de1 = st.visualize_wavefield_2d(output_file, "displacement")
ax.tricontourf(t, de1[-1, :])


# Generate the mesh.
# Elastic moduli and density
m = get_basic_mesh(2)

# Attach parameter to the nodes of each element.
mu = rho * vs**2
kap = rho * (vp**2 - 4 / 3 * vs**2)
par_template = np.ones_like(m.get_element_nodes()[:, :, 0])
m.attach_field("KAPPA", par_template * kap)
m.attach_field("MU", par_template * mu)
m.attach_field("RHO", par_template * rho)
m.attach_field("fluid", np.zeros(m.nelem))

# Attach the mesh and set some custom output.
w.set_mesh(m)
w.output.volume_data.fields = ["displacement"]
w.physics.wave_equation.point_source = [src_vector_2d]

# Run the solver.
output_folder = Path(PROJECT_DIR , "elastic_kappamurho")
output_file = output_folder / "output.h5"
run_salvus(input_file=w, output_folder=output_folder)

# Visualize the results.
f, ax = plt.subplots(1, 1)
ax.set_aspect("equal")
t, de2 = st.visualize_wavefield_2d(output_file, "displacement")
ax.tricontourf(t, de2[-1, :])


# All parameterizations should have produced the same output.
np.testing.assert_allclose(de0, de1, atol=1e-7)
np.testing.assert_allclose(de1, de2, atol=1e-7)


# 3D domains 
w = sn.simple_config.simulation.Waveform()

w.domain.dimension = 3
w.output.volume_data.format = "hdf5"
w.output.volume_data.filename = "output.h5"
w.output.volume_data.sampling_interval_in_time_steps = 100

# Velocities and density
# Generate the mesh.
m = get_basic_mesh(3)

# Attach parameter to the nodes of each element.
par_template = np.ones_like(m.get_element_nodes()[:, :, 0])
m.attach_field("VP", par_template * vp)
m.attach_field("VS", par_template * vs)
m.attach_field("RHO", par_template * rho)
m.attach_field("fluid", np.zeros(m.nelem))

# # Attach the mesh and set some custom output.
w.set_mesh(m)
w.output.volume_data.fields = ["displacement"]
w.physics.wave_equation.point_source = [src_vector_3d]

# # Run the solver.
output_folder = Path(PROJECT_DIR , "elastic_vpvsrho")
output_file = output_folder / "output.h5"
run_salvus(input_file=w, output_folder=output_folder, overwrite=True)

# Visualize the results.
with h5py.File(output_file, "r") as fh:
    de0 = fh["/volume/displacement"][:]
    

# Generate the mesh.
# Lame parameters and density
m = get_basic_mesh(3)

# Attach parameter to the nodes of each element.
mu = rho * vs**2
lam = rho * vp**2 - 2 * mu
par_template = np.ones_like(m.get_element_nodes()[:, :, 0])
m.attach_field("LAMBDA", par_template * lam)
m.attach_field("MU", par_template * mu)
m.attach_field("RHO", par_template * rho)
m.attach_field("fluid", np.zeros(m.nelem))

# # Attach the mesh and set some custom output.
w.set_mesh(m)
w.output.volume_data.fields = ["displacement"]
w.physics.wave_equation.point_source = [src_vector_3d]

# Run the solver.
output_folder = Path("elastic_lambdamurho")
output_file = output_folder / "output.h5"
run_salvus(input_file=w, output_folder=output_folder, overwrite=True)

# # Visualize the results.
with h5py.File(output_file, "r") as fh:
    de1 = fh["/volume/displacement"][:]


# Elastic moduli and density
# Generate the mesh.
m = get_basic_mesh(3)

# Attach parameter to the nodes of each element.
mu = rho * vs**2
kap = rho * (vp**2 - (4 / 3) * vs**2)
par_template = np.ones_like(m.get_element_nodes()[:, :, 0])
m.attach_field("KAPPA", par_template * kap)
m.attach_field("MU", par_template * mu)
m.attach_field("RHO", par_template * rho)
m.attach_field("fluid", np.zeros(m.nelem))

# Attach the mesh and set some custom output.
w.set_mesh(m)
w.output.volume_data.fields = ["displacement"]
w.physics.wave_equation.point_source = [src_vector_3d]

# Run the solver.
output_folder = Path("elastic_kappamurho")
output_file = output_folder / "output.h5"
run_salvus(input_file=w, output_folder=output_folder, overwrite=True)

# Visualize the results.
with h5py.File(output_file, "r") as fh:
    de2 = fh["/volume/displacement"][:]

# All parameterizations should have produced the same output.
np.testing.assert_allclose(de0, de1, atol=1e-7)
np.testing.assert_allclose(de1, de2, atol=1e-7)

