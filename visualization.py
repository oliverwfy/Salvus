import os

SALVUS_FLOW_SITE_NAME = 'oliver_mac'
PROJECT_DIR = "/Users/Oliver/Library/CloudStorage/OneDrive-UniversityofBristol/Salvus_project/project"

import matplotlib.pyplot as plt
import numpy as np
import pathlib
import xarray as xr

import salvus.namespace as sn
from salvus.project.tools.processing import block_processing
from salvus.toolbox.helpers.wavefield_output import (
    WavefieldOutput,
    wavefield_output_to_xarray,
)


print("Opening existing project.")
p = sn.Project(path=PROJECT_DIR)



# run simulation_3 again to get extra outputs.
p.simulations.launch(
    ranks_per_job=4,
    site_name=SALVUS_FLOW_SITE_NAME,
    events=p.events.list(),
    simulation_configuration="simulation_3",
    extra_output_configuration={
        "volume_data": {
            "sampling_interval_in_time_steps": 10,
            "fields": ["displacement", "gradient-of-displacement"],
        },
    },
    # We have previously simulated the same event but without
    # extra output. We have to thus overwrite the existing
    # simulation.
    delete_conflicting_previous_results=True,
)

p.simulations.query(block=True)


p.simulations.get_simulation_output_directory(
    simulation_configuration="simulation_3", event=p.events.list()[0]
)

# use XDMF Reader in Paraview to open the volume_data_output_elastic_volume.xdmf in the directory obtained above.

# get the outputs of receivers 
ed = p.waveforms.get("simulation_3", p.events.list())

# first event, first 

print(
    ed[0].get_receiver_data(
        receiver_name="XX.REC1.", receiver_field="displacement"
    )
)


t_pts = 250

# re-sample with specific sampling rate and number of points
ed[0].set_temporal_interpolation(
    start_time_in_seconds=0.0, sampling_rate_in_hertz=500, npts=t_pts
)


print(
    ed[0].get_receiver_data(
        receiver_name="XX.REC1.", receiver_field="displacement"
    )
)

ed[0].get_receiver_data(
    receiver_name="XX.REC1.", receiver_field="displacement"
).plot()

trace_data = (
    ed[0]
    .get_receiver_data(
        receiver_name="XX.REC1.", receiver_field="displacement"
    )[1]
    .data
)



# times: time step labels 1D array with length t_pts 
# values: outputs of all receivers 2D array (#Rx x t_pts)
times, values = ed[0].get_data_cube(
    receiver_field="displacement", component="Y"
)


plt.plot(times, values[0, :], label="get_data_cube")
plt.plot(times, trace_data, "--", label="get_receiver_data")
plt.legend()
plt.show()


# Salvus uses unstructured meshes,  
# thus post-processing is required to obtain the wavefield output on a regular grid or at a predefined set of points 

# First, retrieve the volumetric output file
output_file = pathlib.Path(
    p.simulations.get_simulation_output_directory("simulation_3", "event_0"),
    "volume_data_output.h5",
)

# Create a WavefieldOutput object for resampling
volume_output = WavefieldOutput.from_file(
    output_file, "displacement", "volume"
)

# output is an xarray.DataArray with with dimensions for time, components and the array of points provided
sampled_output_1 = wavefield_output_to_xarray(
    volume_output, p.simulations.get_mesh("simulation_3").points
).T


x_new = np.linspace(0.0, 2000.0, 200)
y_new = np.linspace(0.0, 1000.0, 100)

# To sample wavefield on a regular grid
sampled_output_2 = wavefield_output_to_xarray(
    volume_output,
    [
        x_new,
        y_new,
    ],
).T


ax = sampled_output_2.sel(c=1, t=sampled_output_2.t[20]).plot(
    shading="gouraud", infer_intervals=False
)
ax.axes.set_aspect("equal")




# resample the output in time 
t_pts_new = 2501
t_range_new = np.linspace(0.0, 0.5, t_pts_new)
field_shape = sampled_output_2.data.shape
new_field_shape = field_shape[0:3] + (t_pts_new,)


resampled_output = block_processing.resample(
    np.concatenate(sampled_output_2.data, axis=0),
    sampled_output_2.t.data,
    t_range_new
    ).reshape(new_field_shape)


idx = np.where(t_range_new == np.round(sampled_output_2.t[20].values, 4))[0]

resampled_xarr = xr.DataArray(
    data = resampled_output.T,
    coords={"t" : t_range_new, "c": [0, 1], "x": x_new, "y": y_new},
    ).T



ax = resampled_xarr.sel(c=1, t=t_range_new[idx]).plot(
    shading="gouraud", infer_intervals=False
)
ax.axes.set_aspect("equal")
