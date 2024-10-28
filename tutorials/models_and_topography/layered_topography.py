import os
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

import salvus.namespace as sn
from salvus.mesh.simple_mesh import rho_from_gardners, vs_from_poisson

SALVUS_FLOW_SITE_NAME = 'oliver_mac'
PROJECT_DIR = "/Users/Oliver/Library/CloudStorage/OneDrive-UniversityofBristol/Salvus_project/model_and_topography"



# x extent of our model 0-5km
x_min = 0.0
x_max = 5000.0

# define a series of 5 discontinuities by specifying their interpolation points in each coordinate
layers_x = [
    np.array([0.0, 0.2, 1.0, 2.0, 3.0, 4.0, 4.8, 5.0]) * 1000,
    np.array([0.0, 0.2, 1.0, 2.0, 3.0, 4.0, 4.8, 5.0]) * 1000,
    np.array([0.0, 0.2, 1.0, 2.0, 3.0, 4.0, 4.8, 5.0]) * 1000,
    np.array([0.0, 1.5, 3.5, 5.0]) * 1000,
    np.array([0.0, 2.5, 5.0]) * 1000,
]

layers_y = [
    np.array([2.0, 2.0, 1.9, 1.7, 2.0, 2.1, 2.0, 2.0]) * 1000,
    np.array([1.6, 1.6, 1.5, 1.4, 1.3, 1.4, 1.5, 1.5]) * 1000,
    np.array([0.5, 0.5, 0.7, 0.6, 1.1, 0.9, 1.2, 1.2]) * 1000,
    np.array([0.2, 0.2, 0.4, 0.4]) * 1000,
    np.array([0.0, 0.0, 0.0]) * 1000,
]

# Define p-velocities.
vp = np.array([2000.0, 2500.0, 2800.0, 3200.0])

# Compute vs and rho.
vs = vs_from_poisson(vp)
rho = rho_from_gardners(vp)

interpolation_styles = [
    "quadratic",
    "quadratic",
    "quadratic",
    "linear",
    "linear",
]

# 4 separated region, 5 interfaces for 2D domain
splines = sn.toolbox.get_interpolating_splines(
    layers_x, layers_y, kind=interpolation_styles
)

# Plot the interfaces.
f = plt.figure(figsize=(10, 5))
x_plot = np.linspace(x_min, x_max)
for top, bot in splines:
    plt.plot(x_plot, top(x_plot))
    plt.plot(x_plot, bot(x_plot))

plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Interfaces")


# Maximum frequency to resolve with elements_per_wavelength.
max_frequency = 20.0

# Generate the mesh
mesh, bnd = sn.toolbox.generate_mesh_from_splines_2d(
    x_min=0,
    x_max=x_max,
    splines=splines,
    elements_per_wavelength=2,
    maximum_frequency=max_frequency,
    use_refinements=True,
    slowest_velocities=vs,
    absorbing_boundaries=(["x0", "x1", "y0"], 10.0),
)

# sum the mesh regions into one continuous mesh
mesh = np.sum(mesh)

# Add info about absorbing boundaries
mesh.attach_global_variable("max_dist_ABC", bnd)
mesh.attach_global_variable("ABC_side_sets", ", ".join(["x0", "x1", "y0"]))
mesh.attach_global_variable("ABC_vel", min(vs))
mesh.attach_global_variable("ABC_freq", max_frequency / 2.0)
mesh.attach_global_variable("ABC_nwave", 5.0)

mesh  # Visualize the mesh.


nodes = mesh.get_element_nodes()[:, :, 0]
vp_a, vs_a, ro_a = np.ones((3, *nodes.shape))
for _i, (vp_val, vs_val, ro_val) in enumerate(zip(vp, vs, rho)):
    # Find which elements are in a given region.
    idx = np.where(mesh.elemental_fields["region"] == _i)

    # Set parameters in that region to a constant value.
    vp_a[idx] = vp_val
    vs_a[idx] = vs_val
    ro_a[idx] = ro_val

# Attach parameters.
for k, v in zip(["VP", "VS", "RHO"], [vp_a, vs_a, ro_a]):
    mesh.attach_field(k, v)

# Attach acoustic / elastic flag.
mesh = sn.toolbox.detect_fluid(mesh)

mesh


p = sn.Project.from_domain(
    path=Path(PROJECT_DIR , "Project_layered_topography"),
    domain=sn.domain.dim2.BoxDomain(x0=x_min, x1=x_max, y0=0.0, y1=2700.0),
    load_if_exists=True,
)
p.viz.nb.domain()

# Location of source (this will snap to the closet side-set).
loc = [x_max / 2, 2500.0]

# Create the source.
source = sn.simple_config.source.cartesian.SideSetVectorPoint2D(
    fx=1,
    fy=1,
    point=loc,
    direction="y",
    side_set_name="y1",
)

receivers = (
    sn.simple_config.receiver.cartesian.SideSetHorizontalPointCollection2D(
        x=np.linspace(x_min, x_max, 1000),
        offset=-1.0,
        side_set_name="y1",
        station_code="xx",
        fields=["velocity", "strain"],
    )
)

p.add_to_project(
    sn.Event(event_name="event", sources=source, receivers=receivers)
)


if "sim" not in p.simulations.list():
    p.add_to_project(
        sn.UnstructuredMeshSimulationConfiguration(
            name="sim",
            unstructured_mesh=mesh,
            event_configuration=sn.EventConfiguration(
                wavelet=sn.simple_config.stf.Ricker(
                    center_frequency=max_frequency / 2
                ),
                waveform_simulation_configuration=sn.WaveformSimulationConfiguration(
                    end_time_in_seconds=2.0
                ),
            ),
        ),
    )


p.viz.nb.simulation_setup("sim", events=["event"])

p.simulations.launch(
    simulation_configuration="sim",
    events="event",
    site_name=SALVUS_FLOW_SITE_NAME,
    ranks_per_job=4,
)
p.simulations.query(block=True)

waveform_data = p.waveforms.get(data_name="sim", events="event")[0]

waveform_data.plot(
    component="X",
    receiver_field="velocity",
    plot_types=["shotgather"],
)
waveform_data.plot(
    component="Y",
    receiver_field="velocity",
    plot_types=["shotgather"],
)

waveform_data.plot(
    component="0",
    receiver_field="strain",
    plot_types=["shotgather"],
)
waveform_data.plot(
    component="1",
    receiver_field="strain",
    plot_types=["shotgather"],
)
waveform_data.plot(
    component="2",
    receiver_field="strain",
    plot_types=["shotgather"],
)