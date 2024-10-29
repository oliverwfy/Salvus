# Standard Python packages
import os
from pathlib import Path

# Third-party imports.
import numpy as np
import matplotlib.pyplot as plt
import pyasdf

# Salvus packages
import salvus.flow.api
import salvus.flow.simple_config as sc
from salvus.mesh.simple_mesh import basic_mesh

SALVUS_FLOW_SITE_NAME = 'oliver_mac'
PROJECT_DIR = "/Users/Oliver/Library/CloudStorage/OneDrive-UniversityofBristol/Salvus_project/sources_and_receivers"

# A single array. (source time function)
stf = sc.stf.Custom.from_array(
    np.sin(np.linspace(0, 4 * np.pi, 100)),
    sampling_rate_in_hertz=2.0,
    start_time_in_seconds=0.0,
)

stf.plot()

# Combine with a source object to create a complete source object.
# Note the different weights (determined by fx and fy) will define the final source together
# with the STF.
src = sc.source.cartesian.VectorPoint2D(
    x=10.0, y=0.0, fx=1e-2, fy=-2e-2, source_time_function=stf
)


# It is also possible to specify a separate array for every components.
array = np.array(
    [
        1.5 * np.sin(np.linspace(0, 4 * np.pi, 100)),
        -3.0 * np.sin(np.linspace(0, 4 * np.pi, 100)),
    ]
)

stf = sc.stf.Custom.from_array(
    array, sampling_rate_in_hertz=2.0, start_time_in_seconds=-10.0
)
stf.plot()

# Note that in this case the weights should be set to 1.0
src = sc.source.cartesian.VectorPoint2D(
    x=10.0, y=0.0, fx=1.0, fy=1.0, source_time_function=stf
)



# Ricker wavelet
wavelet_width_in_seconds = 0.1075
time_step_in_seconds = 1e-3
center_frequency = 14.5

sigma_2 = 1 / (np.pi * center_frequency) ** 2

time = np.linspace(
    -wavelet_width_in_seconds,
    wavelet_width_in_seconds,
    int((2 * wavelet_width_in_seconds / time_step_in_seconds)),
)

sampling_rate_in_hertz = 1.0 / time_step_in_seconds

wavelet = (1 - (2 * time**2) / sigma_2) * np.exp(-(time**2) / sigma_2)

# plot the wavelet
plt.plot(time, wavelet)
plt.xlabel("time [s]")
plt.ylabel("amplitude")
plt.title("Custom ricker wavelet")
plt.show()


# Simple 2D elastic mesh.
mesh = basic_mesh.CartesianHomogeneousIsotropicElastic2D(
    vp=3200.0,
    vs=1847.5,
    rho=2200.0,
    x_max=2000.0,
    y_max=1000.0,
    max_frequency=25.0,
).create_mesh()


receiver = sc.receiver.cartesian.Point2D(
    x=1400.0, y=700.0, station_code="XX", fields=["displacement"]
)


# Spatial weights of the vectorial source.
fx = 1e-5
fy = -0.8e-4

# Location.
sx = 1000.0
sy = 500.0

# Option 1 - Parameterized STF.
stf_1 = sc.stf.Ricker(center_frequency=14.5)
source_1 = custom_source = sc.source.cartesian.VectorPoint2D(
    x=sx, y=sy, fx=fx, fy=fy, source_time_function=stf_1
)

# Option 2 - single-component STF and associated weights.
stf_2 = sc.stf.Custom.from_array(
    array=wavelet,
    sampling_rate_in_hertz=sampling_rate_in_hertz,
    start_time_in_seconds=time[0],
)
source_2 = sc.source.cartesian.VectorPoint2D(
    x=sx, y=sy, fx=fx, fy=fy, source_time_function=stf_2
)

# Option 3 - multi-component STF and unit weights.
source_time_function = [wavelet * fx, wavelet * fy]
stf_3 = sc.stf.Custom.from_array(
    array=source_time_function,
    sampling_rate_in_hertz=sampling_rate_in_hertz,
    start_time_in_seconds=time[0],
)
source_3 = sc.source.cartesian.VectorPoint2D(
    x=sx, y=sy, fx=1.0, fy=1.0, source_time_function=stf_3
)


stf_1.plot()
stf_2.plot()
stf_3.plot()



for _i, src in enumerate([source_1, source_2, source_3]):
    sim = sc.simulation.Waveform(mesh=mesh, sources=src, receivers=receiver)
    sim.physics.wave_equation.end_time_in_seconds = 0.5

    salvus.flow.api.run(
        site_name=SALVUS_FLOW_SITE_NAME,
        input_file=sim,
        ranks=4,
        output_folder=Path(PROJECT_DIR, f"output_custom_{_i}"),
        overwrite=True,
    )
    
_, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))

for _i in range(3):
    folder = Path(PROJECT_DIR, f"output_custom_{_i}")
    with pyasdf.ASDFDataSet(folder / "receivers.h5") as ds:
        for _j, ax in enumerate(axes):
            tr = ds.waveforms.XX_XX.displacement[_j]
            ax.plot(tr.times(), tr.data, label=f"Source {_i}")
            ax.set_title(f"Component {_j + 1}")
axes[0].legend()
axes[1].legend()
plt.show()