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
ranks = 4


mesh = sn.simple_mesh.basic_mesh.CartesianHomogeneousAcoustic2D(
    vp=1000.0, rho=1000.0, x_max=1.0, y_max=1.0, max_frequency=12500.0
).create_mesh()



mesh.elemental_fields.clear()

mesh_order_1 = mesh.copy()
mesh_order_1.change_tensor_order(1)

mesh_order_4 = mesh.copy()
mesh_order_4.change_tensor_order(4)

# axis = np.random.choice([0, 1])
axis = 1
coords_order_1 = mesh_order_1.get_element_nodes()[:, :, axis].flatten()
coords_order_4 = mesh_order_4.get_element_nodes()[:, :, axis].flatten()

def li(x: np.ndarray, min_val, max_val) -> np.ndarray:
    """Linearly interpolate a parameter over a coordinate array.

    Parameters
    ----------
    x : np.ndarray
        Flattened coordinate array (either x or y).
    min_val : float
        Value of the parameter at the minimum coordinate.
    max_val : float
        Value of the parameter at the maximum coordinate.

    Returns
    -------
    np.ndarray
        An array of values linearly interpolated over the coordinate range.

    """
    m = (max_val - min_val) / (np.max(x) - np.min(x))
    return m * x + min_val


vp0, vp1 = 1000, 2000
rho = 1000
m10 = m11 = 1 / rho
m00 = m10 / vp0**2
m01 = m11 / vp1**2


# attach the desired parameter to the 1st order mesh
m0 = (np.apply_along_axis(li, 0, coords_order_1, m00, m01)).reshape(
    mesh_order_1.nelem, mesh_order_1.nodes_per_element
) * 0.5

m1 = (np.apply_along_axis(li, 0, coords_order_1, m10, m11)).reshape(
    mesh_order_1.nelem, mesh_order_1.nodes_per_element
) * 0.5

mesh_order_1.attach_field("M0", m0)
mesh_order_1.attach_field("M1", m1)
mesh_order_1.attach_field("fluid", np.ones(mesh.nelem))
mesh_order_1

# attach the desired parameter to the 4th order mesh
m0 = (np.apply_along_axis(li, 0, coords_order_4, m00, m01)).reshape(
    mesh_order_4.nelem, mesh_order_4.nodes_per_element
) * 0.5

m1 = (np.apply_along_axis(li, 0, coords_order_4, m10, m11)).reshape(
    mesh_order_4.nelem, mesh_order_4.nodes_per_element
) * 0.5

mesh_order_4.attach_field("M0", m0)
mesh_order_4.attach_field("M1", m1)
mesh_order_4.attach_field("fluid", np.ones(mesh.nelem))

mesh_order_4

# set up simulations
# Source in the left borehole.
src = sn.simple_config.source.cartesian.ScalarPoint2D(
    x=0.1,
    y=0.5,
    f=1,
    source_time_function=sn.simple_config.stf.Ricker(center_frequency=1e4),
)

# String of receivers in the right borehole.
recs = sn.simple_config.receiver.cartesian.SideSetVerticalPointCollection2D(
    y=np.linspace(0, 1, 101),
    side_set_name="x1",
    station_code="xx",
    fields=["phi"],
    offset=-0.1,
)

w1 = sn.simple_config.simulation.Waveform(
    mesh=mesh_order_1, sources=src, receivers=recs
)


# Output receivers in simplified HDF5 format.
w1.output.point_data.format = "hdf5"

# Run the simulation long enough that wave hit the receivers.
w1.physics.wave_equation.end_time_in_seconds = 1e-3

# Output the state variable at all time steps (for comparison).
w1.output.volume_data.format = "hdf5"
w1.output.volume_data.fields = ["phi"]
w1.output.volume_data.filename = "output.h5"
w1.output.volume_data.sampling_interval_in_time_steps = 1

# Validate and plot.
w1.validate()
# w1


of1 = PROJECT_DIR / Path("output_order_1")

# Run.
sn.api.run(
    ranks=ranks,
    get_all=True,
    input_file=w1,
    overwrite=True,
    site_name=SALVUS_FLOW_SITE_NAME,
    output_folder=of1,
)

ed = sn.EventData.from_output_folder(of1)
ed.plot(receiver_field="phi", component="A")


# Get output data into notebook.
f, ax = plt.subplots(1, 2, figsize=(15, 5))
tri, d1_vol = st.visualize_wavefield_2d(of1 / "output.h5", field="phi")
d1, t_step, extent = st.get_shotgather(of1 / "receivers.h5", field="phi", axis=1)

# Plot shotgather.
ax[0].imshow(d1, extent=extent, aspect="auto")
ax[0].set_xlabel("Position of receiver in y axis (m)")
ax[0].set_ylabel("Time (s)")
ax[0].set_title("Shotgather")

# Plot last time step of wavefield.
ax[1].tricontourf(tri, d1_vol[100, :])
ax[1].set_xlabel(" x axis (m)")
ax[1].set_ylabel(" y axis (m)")
ax[1].set_title("Wavefield snapshot")




# set up order 4 mesh
w4 = sn.simple_config.simulation.Waveform(
    mesh=mesh_order_4, sources=src, receivers=recs
)

w4.output.point_data.format = "hdf5"
w4.physics.wave_equation.end_time_in_seconds = 1e-3

w4.output.volume_data.format = "hdf5"
w4.output.volume_data.fields = ["phi"]
w4.output.volume_data.filename = "output.h5"
w4.output.volume_data.sampling_interval_in_time_steps = 1

w4.validate()

# Output folder 4.
of4 = PROJECT_DIR / Path("output_order_4")

# Run.
sn.api.run(
    input_file=w4,
    site_name=SALVUS_FLOW_SITE_NAME,
    ranks=ranks,
    output_folder=of4,
    overwrite=True,
    get_all=True,
)

f = plt.figure(figsize=(15, 10))
gs = f.add_gridspec(2, 2)
ax0 = f.add_subplot(gs[0, 0])
ax1 = f.add_subplot(gs[0, 1])
ax2 = f.add_subplot(gs[1, :])

d1, dt, extent = st.get_shotgather(of1 / "receivers.h5", field="phi", axis=1)
d4, dt, extent = st.get_shotgather(of4 / "receivers.h5", field="phi", axis=1)


ax0.imshow(d1, extent=extent, aspect="auto")
ax0.set_xlabel("Position (m)")
ax0.set_ylabel("Time (s)")
ax0.set_title("Shotgather (order 1)")

ax1.imshow(d4, extent=extent, aspect="auto")
ax1.set_xlabel("Position (m)")
ax1.set_ylabel("Time (s)")
ax1.set_title("Shotgather (order 4)")


ax2.set_title("Trace comparison of a specific receiver (#50)")
ax2.plot(np.arange(d4.shape[0]) * dt, d1[:, 50], label="Order 1")
ax2.plot(np.arange(d4.shape[0]) * dt, d4[:, 50], label="Order 4", ls="dashed")
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Amplitude")
ax2.legend()


np.testing.assert_allclose(d1, d4, atol=1e-5 * np.max(np.abs(d1)))
