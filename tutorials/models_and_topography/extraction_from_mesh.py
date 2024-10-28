import numpy as np
import xarray as xr
import salvus.namespace as sn
from salvus.mesh.unstructured_mesh_utils import extract_model_to_regular_grid

from pathlib import Path


SALVUS_FLOW_SITE_NAME = 'oliver_mac'
PROJECT_DIR = "/Users/Oliver/Library/CloudStorage/OneDrive-UniversityofBristol/Salvus_project/model_and_topography"

# Gaussian
x, y = np.linspace(0.0, 1.0, 101), np.linspace(0.0, 1.0, 101)
xx, yy = np.meshgrid(x, y, indexing="xy")
g = np.exp(-(((xx - 0.5) ** 2 + (yy - 0.5) ** 2) / (2 * 0.2**2)))

# Material parameters
vp = 2 * g + 1
rho = vp / 2

# Xarray dataset
ds_model_2d = xr.Dataset(
    coords={"x": x, "y": y},
    data_vars={"VP": (["x", "y"], vp), "RHO": (["x", "y"], rho)},
)

# Salvus wrapper.
m = sn.model.volume.cartesian.GenericModel(name="blob", data=ds_model_2d)

# Plot
m.ds.VP.plot()

p_2d = sn.Project.from_volume_model(Path(PROJECT_DIR,"Proj_slice_2d"), m, True)

sc = sn.SimulationConfiguration(
    name="blob_fmax",
    event_configuration=None,
    elements_per_wavelength=2.0,
    max_frequency_in_hertz=10.0,
    model_configuration=sn.ModelConfiguration(
        background_model=None, volume_models="blob"
    ),
)

p_2d.add_to_project(sc, overwrite=True)

mesh_0 = p_2d.simulations.get_mesh("blob_fmax")
mesh_0


# Extract to regular grid from mesh
# define coordinates 
ds_extract_2d = xr.Dataset(
    coords={"x": np.linspace(0, 1.0, 101), "y": np.linspace(0, 1.0, 101)}
)

ds_extract_2d = extract_model_to_regular_grid(
    mesh_0, ds_extract_2d, ["VP", "RHO"]
)

ds_extract_2d.VP.plot()


# Extract profile at y = 0.5
ds_extract_2d_line = xr.Dataset(
    coords={"x": np.linspace(0.0, 1.0, 101), "y": [0.5]}
)
ds_extract_2d_line = extract_model_to_regular_grid(
    mesh_0, ds_extract_2d_line, ["VP", "RHO"]
)

ds_extract_2d_line.VP.plot()


# 3D 
# Gaussian
x, y, z = (
    np.linspace(0.0, 1.0, 101),
    np.linspace(0.0, 1.0, 101),
    np.linspace(0.0, 1.0, 101),
)

xx, yy, zz = np.meshgrid(x, y, z, indexing="xy")
g = np.exp(
    -(((xx - 0.5) ** 2 + (yy - 0.5) ** 2 + (zz - 0.5) ** 2) / (2 * 0.2**2))
)

# Pars
vp_3d = 2 * g + 1
rho_3d = vp_3d / 2

# Xarray dataset
ds_model_3d = xr.Dataset(
    coords={"x": x, "y": y, "z": z},
    data_vars={
        "VP": (["x", "y", "z"], vp_3d),
        "RHO": (["x", "y", "z"], rho_3d),
    },
)

# Salvus wrapper.
m_3d = sn.model.volume.cartesian.GenericModel(name="blob", data=ds_model_3d)

p_3d = sn.Project.from_volume_model(Path(PROJECT_DIR,"Proj_slice_3d"), m_3d, True)

p_3d.add_to_project(sc, overwrite=True)

mesh_3d = p_3d.simulations.get_mesh("blob_fmax")
mesh_3d

ds_extract_3d = xr.Dataset(
    coords={
        "x": np.linspace(0, 1.0, 101),
        "y": np.linspace(0, 1.0, 101),
        "z": np.linspace(0, 1.0, 11),
    }
)

ds_extract_3d = extract_model_to_regular_grid(
    mesh_3d, ds_extract_3d, ["VP", "RHO"], verbose=True
)

# Get an xy slice at z == 0.5.
ds_extract_3d.VP.sel(z=0.5).plot()