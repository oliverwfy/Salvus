import os
from salvus.flow import simple_config
from pathlib import Path
import salvus.flow.api
from salvus.mesh import simple_mesh


SALVUS_FLOW_SITE_NAME = 'oliver_mac'
PROJECT_DIR = "/Users/Oliver/Library/CloudStorage/OneDrive-UniversityofBristol/Salvus_project/sources_and_receivers"


# receivers
rec = simple_config.receiver.cartesian.Point2D(
    # Cartesian coordinates.
    x=2000.0,
    y=2000.0,
    # The network is optional but helps to group receivers.
    network_code="XX",
    # The name of the receiver.
    station_code="A1",
    # An additional level to group receivers.
    location_code="",
    # At least one output field is required. More are possible.
    # Have a look at the API documentation for a list of all
    # available fields.
    fields=["displacement", "velocity", "acceleration"],
)


# They are internally represented as dictionaries exactly
# corresponding to what SalvusCompute demands.
print(rec)


src = simple_config.source.cartesian.VectorPoint2D(
    # Coordinates of the source.
    x=500.0,
    y=1000.0,
    # Force vector in x and y direction in N.
    fx=1e5,
    fy=-1e4,
    # It also requires a source time function.
    source_time_function=simple_config.stf.Ricker(center_frequency=1.0),
)

# They are again internally represented as a dictionary.
print(src)



m = simple_mesh.CartesianHomogeneousIsotropicElastic2D(
    vp=2000.0,
    vs=1500.0,
    rho=2000.0,
    x_max=3000.0,
    y_max=2000.0,
    max_frequency=2.0,
)

w = simple_config.simulation.Waveform(mesh=m.create_mesh())
w.add_receivers(rec)
w.add_sources(src)
w


salvus.flow.api.run(
    site_name=SALVUS_FLOW_SITE_NAME, input_file=w, output_folder=Path(PROJECT_DIR ,"output"), 
    overwrite=True
)