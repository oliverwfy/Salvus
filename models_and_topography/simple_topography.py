import os
import numpy as np
import xarray as xr
import salvus.namespace as sn
from pathlib import Path


SALVUS_FLOW_SITE_NAME = 'oliver_mac'
PROJECT_DIR = "/Users/Oliver/Library/CloudStorage/OneDrive-UniversityofBristol/Salvus_project/model_and_topography"


d = sn.domain.dim2.BoxDomain(x0=0, x1=20e3, y0=-5e3, y1=0)
d.plot()
p = sn.Project.from_domain(path=Path(PROJECT_DIR, "Project_simple_topography"), 
                           domain=d, load_if_exists=True)

max_x = d.bounding_box[:, 0].max()
x = np.linspace(0, max_x, 10000)
y = np.cos(x * 2 * np.pi / max_x) * 1e3 


ds = xr.Dataset(
    data_vars={"dem": (["x"], y, {"reference_elevation": 0.0})},
    coords={"x": x},
)

tm = sn.topography.cartesian.SurfaceTopography(name="surface", data=ds)
tm.ds.dem.plot()

p.add_to_project(tm, overwrite=True)

vp, vs, rho = 4e3, 2.4e3, 2.6e3
bm = sn.model.background.homogeneous.IsotropicElastic(vp=vp, vs=vs, rho=rho)

src = sn.simple_config.source.cartesian.SideSetVectorPoint2D(
    point=(5.0e3, 0.0),
    direction="y",
    offset=-10.0,
    side_set_name="y1",
    fx=0.0,
    fy=1.0,
)

recs = [
    sn.simple_config.receiver.cartesian.SideSetPoint2D(
        point=(x, 0.0),
        direction="y",
        offset=0.0,
        side_set_name="y1",
        station_code=f"XX_{i}",
        fields=["velocity", "strain"],
    )
    for i, x in enumerate(np.linspace(12.5e3, 17.5e3, 11))
]

p.add_to_project(sn.EventCollection.from_sources(sources=src, receivers=recs))

f_max = 10.0
p.add_to_project(
    sn.SimulationConfiguration(
        tensor_order=2,
        name="sim_with_topo",
        elements_per_wavelength=2.0,
        max_frequency_in_hertz=f_max,
        model_configuration=sn.ModelConfiguration(background_model=bm),
        topography_configuration=sn.TopographyConfiguration("surface"),
        event_configuration=sn.EventConfiguration(
            sn.simple_config.stf.Ricker(center_frequency=f_max / 2),
            sn.WaveformSimulationConfiguration(end_time_in_seconds=10.0),
        ),
        absorbing_boundaries=sn.AbsorbingBoundaryParameters(
            reference_velocity=vp,
            number_of_wavelengths=0.0,
            reference_frequency=f_max / 2,
        ),
    ),
    overwrite=True,
)


p.viz.nb.simulation_setup("sim_with_topo", events="event_0000")


p.simulations.launch(
    "sim_with_topo",
    events="event_0000",
    site_name=SALVUS_FLOW_SITE_NAME,
    ranks_per_job=4,
)

p.simulations.query(block=True)

p.waveforms.get(
    data_name="SYNTHETIC_DATA:sim_with_topo", events=["event_0000"]
)[0].plot(component="X", receiver_field="velocity")
