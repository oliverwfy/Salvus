from material_ela_constants.Elastic_Material import Austenite
import salvus.namespace as sn
from salvus.material import elastic
from pathlib import Path
import salvus
import numpy as np
from my_code.utilities import *
from salvus.flow.simple_config.receiver.cartesian import Point2D, Point3D
from salvus.flow.simple_config.source.cartesian import VectorPoint2D, VectorPoint3D
import salvus.mesh.layered_meshing as lm
from datetime import datetime
import matplotlib.pyplot as plt

# Salvus site name
SALVUS_FLOW_SITE_NAME = 'oliver_wsl'
RANKS_PER_JOB = 4

# Wokring dir
WORKING_DIR = '/home/oliver/workspace/Salvus/elastic_model/anisotropic/'

# Directories in WSL
PROJECT_DIR = '/home/oliver/workspace/Salvus/elastic_model/anisotropic/Project'
IMAGE_DIR = '/home/oliver/workspace/Salvus/elastic_model/anisotropic/image'
DATA_DIR = '/home/oliver/workspace/Salvus/elastic_model/anisotropic/data'

# Directories in Windows
PROJECT_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/Project'
DATA_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/data'
IMAGE_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/image'


# create dir if it does not exist
Path(IMAGE_DIR).mkdir(parents=True, exist_ok=True)
Path(DATA_DIR).mkdir(parents=True, exist_ok=True)
Path(IMAGE_DIR_WIN).mkdir(parents=True, exist_ok=True)
Path(DATA_DIR_WIN).mkdir(parents=True, exist_ok=True)



project_name = 'attenuation'



matl = Austenite()      # load material's elasticity tensor
f_c = 2*1e6             # centre freuency     




# 2D box domain parameters (length in m)
x_length = 10 * 1e-3
y_length = 2 * 1e-3
z_length = 10 * 1e-3


x_range = (0., x_length) 
y_range = (0., y_length) 
z_range = (0., z_length) 

# define 2D box domain
domain = sn.domain.dim3.BoxDomain(x0=x_range[0], x1=x_range[1], y0=y_range[0], y1=y_range[1], z0=z_range[0], z1=z_range[1])


# create Project from domain
p = sn.Project.from_domain(path=Path(PROJECT_DIR_WIN, project_name), 
                           domain=domain, load_if_exists=True)


# # number of point vector sources and receivers.
# n_srcs = 64
# n_rxs = 64


# # range of sources and receivers
# srcs_range = np.array([0, 1]) * x_length
# rxs_range = np.array([0.3, 0.7]) * x_length


# # positions of srcs and rxs
# srcs_pos = [Vector(np.round(x, 5), y_range[1]) 
#             for x in np.linspace(*srcs_range, n_srcs)
#             ]
# rxs_pos = [Vector(np.round(x, 5), y_range[0]) 
#            for x in np.linspace(*rxs_range, n_rxs)
#            ]








n_srcs = 1001

srcs_pos = [Vector(x, y_range[1]/2, z_range[1]) for x in np.linspace(0, x_length, n_srcs).]


# rxs_pos = [Vector(x_length*0.2, y_range[1]),   Vector(x_length/2, y_range[1]),   Vector(x_length*0.8, y_range[1]),
#            Vector(x_length*0.2, y_range[1]/2), Vector(x_length/2, y_range[1]/2), Vector(x_length*0.8, y_range[1]/2),
#            Vector(x_length*0.2, y_range[0]), Vector(x_length/2, y_range[0]), Vector(x_length*0.8, y_range[0]),]



src_dirw = Vector(0, 2, 0)      # weights applied in x, y, z directions respevtively. 
fileds = ["displacement"]       # received fileds



# vector source 3D with weights fx and fy in x and y directions, respectively.
srcs = [VectorPoint3D(x=s.x,y=s.y,z=s.z, fx=src_dirw.x, fy=src_dirw.y, fz=src_dirw.z) 
        for s in srcs_pos
        ]




# create events for simulation
events = []

# # (array) add all receivers to each event of one point source
# for i, src in enumerate(srcs):
#     rxs = [Point2D(x=r.x, y=r.y, 
#             station_code=f"REC{i + 1}",
#             # Note that one can specify the desired recording field here.
#             fields=fileds,)
#         for i, r in enumerate(rxs_pos)
#         ]

#     events.append(
#         sn.Event(event_name=f"event_{i}", sources=src, receivers=rxs)
#     )


n_rxs = 101
rxs_pos = [Vector(x_range[1]/2, y_range[1]/2, z) for z in np.linspace(z_length,0, n_rxs)]

# (array) add all receivers to each event of one point source
rxs = [Point3D(x=r.x, y=r.y, z=r.z,
        station_code=f"REC_{i + 1}",
        # Note that one can specify the desired recording field here.
        fields=fileds,)
    for i, r in enumerate(rxs_pos)
    ]


events.append(sn.Event(event_name=f"event_0", sources=srcs, receivers=rxs))


# add the events to Project
add_events_to_Project(p, events)


# material_unoriented = sn.material.elastic.triclinic.TensorComponents.from_params(**matl.rotated_parameters(0))

# rotate 45 degree.
layer_1 = sn.material.elastic.triclinic.TensorComponents.from_params(**matl.rotated_parameters(0))
layer_2 = sn.material.elastic.triclinic.TensorComponents.from_params(**matl.rotated_parameters(np.pi/3))
interface = sn.layered_meshing.interface.Hyperplane.at(z_range[1]/2)


# layered_model = sn.layered_meshing.LayeredModel([layer_1])

# simulation_name = "mesh_unoriented"



layered_model = sn.layered_meshing.LayeredModel([
    layer_1,
    interface,
    layer_2])


simulation_name = "mesh_interface_60"


"""
Absorbing Boundary (free-surface)
"""

reference_velocity = 3000           # wave velocity in the absorbing boundary layer
number_of_wavelengths=1           # number of wavelengths to pad the domain by
reference_frequency = f_c           # reference frequency for the distance calculation
free_surfaces = ['z1']       # free surfaces for absorbing boundary layer


layered_model_abs = sn.layered_meshing.MeshingProtocol(
    layered_model,
    ab=salvus.mesh.simple_mesh.basic_mesh.AbsorbingBoundaryParameters(
        free_surface=free_surfaces,
        number_of_wavelengths=number_of_wavelengths,
        reference_velocity=reference_velocity,
        reference_frequency=reference_frequency,
    ),
)


# Define the mesh resolution using salvus
mesh_resolution = sn.MeshResolution(
    reference_frequency=f_c * 2,  # Reference frequency for the mesh
    elements_per_wavelength= 3.5,  # Number of elements per wavelength
    model_order= 2  # Model order for the mesh
)


# Create the mesh for the grains model using the defined domain and mesh resolution
mesh = sn.layered_meshing.mesh_from_domain(
    domain=domain,
    model=layered_model_abs,
    mesh_resolution=mesh_resolution
)





# layered_model_abs_rotated = sn.layered_meshing.MeshingProtocol(
#     layered_model_rotated,
#     ab=salvus.mesh.simple_mesh.basic_mesh.AbsorbingBoundaryParameters(
#         free_surface=free_surfaces,
#         number_of_wavelengths=number_of_wavelengths,
#         reference_velocity=reference_velocity,
#         reference_frequency=reference_frequency,
#     ),
# )


# # Create the mesh for the grains model using the defined domain and mesh resolution
# mesh_rotated = sn.layered_meshing.mesh_from_domain(
#     domain=domain,
#     model=layered_model_abs_rotated,
#     mesh_resolution=mesh_resolution
# )





# # shape (m x n x d)  m=#ele, n=(tensor_order+1)^dim d = dim
# mess_nodes = mesh.get_element_nodes()       
# mesh.get_element_nodes().shape


# """
# Temporal configuration:
#     source time function
#     simulation time duration
# """

# wavelet (input source time function) 
wavelet = sn.simple_config.stf.Ricker(center_frequency=f_c)


# waveform simulation temporal parameters
start_time = -0.3*1e-6
end_time = 8*1e-6

# waveform simulation configuration
waveform_config = sn.WaveformSimulationConfiguration(
        start_time_in_seconds=start_time,
        end_time_in_seconds=end_time,
        )

# event configuaration
event_config = sn.EventConfiguration(
    wavelet=wavelet,
    waveform_simulation_configuration=waveform_config,
    )



# # save figure of event config 
# fig = event_config.wavelet.plot(show=False)
# plt.savefig(Path(IMAGE_DIR_WIN, 'Ricker.png'))









sim_config = sn.UnstructuredMeshSimulationConfiguration(
    unstructured_mesh=mesh,
    name=simulation_name,
    event_configuration=event_config,
    )



# add simulation configuration to Project
p.add_to_project(
    sim_config, overwrite=True
    )



# visualization of mesh and simulation set-up
p.viz.nb.simulation_setup(
    simulation_configuration=simulation_name, events=p.events.list()
    )


dofs =  mesh.number_of_nodes
print(f'Start simulation: {simulation_name}')
print(f'Dofs (number of nodes): {dofs}')


start_time = datetime.now()


# launch simulations
p.simulations.launch(
    ranks_per_job=RANKS_PER_JOB,
    site_name=SALVUS_FLOW_SITE_NAME,
    events=p.events.list(),
    simulation_configuration=simulation_name,
    delete_conflicting_previous_results=True,
    )




# # simulation with volume data (full wavefield)
# p.simulations.launch(
#     ranks_per_job=RANKS_PER_JOB,
#     site_name=SALVUS_FLOW_SITE_NAME,
#     events=p.events.list(),
#     simulation_configuration=simulation_name,
#     extra_output_configuration={
#         "volume_data": {
#             "sampling_interval_in_time_steps": 20,
#             "fields": ["displacement"],
#         },
#     },
#     # We have previously simulated the same event but without
#     # extra output. We have to thus overwrite the existing
#     # simulation.
#     delete_conflicting_previous_results=True,
# )

p.simulations.query(block=True)


end_time = datetime.now()


execution_time_seconds = (end_time - start_time).total_seconds()
minutes = int(execution_time_seconds // 60)  # Extract minutes
seconds = execution_time_seconds % 60  # Extract remaining seconds

print(f"Execution time: {minutes} minutes and {seconds:.2f} seconds")