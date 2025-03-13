import re
import numpy as np
from pathlib import Path
from collections import namedtuple





def add_events_to_Project(Project, events):
    for event in Project.events.list():
        Project.events.delete(event) 
    for event in events:
        Project.add_to_project(event)
    return None

def reorder_list(lst, order):
    return [lst[i] for i in order]


def extract_numbers(lst):
    return [int(num) for s in lst for num in re.findall(r'\d+', s)]


def reorder_events_list(lst):

    e_nums = extract_numbers(lst)
    order = []
    for i in range(len(e_nums)): 
        order.append(int(np.where(np.array(e_nums) == i)[0][0]))
    return reorder_list(lst, order)

    
def source_location(event_data):
    return [(src.location[0], src.location[1]) for e in event_data for src in e.sources] 

def receriver_location(event_data):
    return [(rx.location[0], rx.location[1]) for rx in event_data[0].receivers]


def time_from_ed(event_data, temporal_interpolation=False):
    
    if temporal_interpolation:
        time = None
    else:
        time = event_data[0].get_waveform_data_xarray('displacement').time.values        
    return time


def fmc_data_from_ed(event_data, save_dir=None):
    time = time_from_ed(event_data)
    time = time.reshape(len(time), -1)
    srcs_loc = source_location(event_data)
    rxs_loc = receriver_location(event_data)
    fmc_data = np.zeros([len(time), 2,len(srcs_loc)*len(rxs_loc)])

    for i, _ in enumerate(event_data):
        fmc_data[:, 0,i*len(rxs_loc):(i+1)*len(rxs_loc)] = \
        event_data[0].get_data_cube(receiver_field='displacement', component='X')[1].T

        fmc_data[:, 1,i*len(rxs_loc):(i+1)*len(rxs_loc)] = \
        event_data[0].get_data_cube(receiver_field='displacement', component='Y')[1].T

    if save_dir:
        np.save(Path(save_dir, "time.npy"), time)
        np.save(Path(save_dir, "fmc_data.npy"), fmc_data)   
        np.save(Path(save_dir, "rxs_loc.npy"), rxs_loc)
        np.save(Path(save_dir, "srcs_loc.npy"), srcs_loc)
        
    return (fmc_data, time, rxs_loc, srcs_loc) 


class Vector:
    def __init__(self, x=None, y=None, z=None):
        self.x = x
        self.y = y
        self.z = z
    
    def distance(self, other):
        """Calculates the Euclidean distance to another Vector."""
        return np.linalg.norm(self.array() - other.array(), ord=2)

    def array(self):
        """Returns the vector as a NumPy array."""
        if self.z: 
            return np.array([self.x, self.y, self.z])
        else:
            return np.array([self.x, self.y])



def compute_slowness(C, rho, theta):
    # Define propagation direction (sin(theta), 0, cos(theta))
    n3 = np.cos(np.radians(theta))
    n1 = np.sin(np.radians(theta))
    # VTI elasticity tensor components
    C11, C12, C13, C33, C44, C66 = C

    # Construct the Christoffel matrix
    Gamma = np.array([
        [C11 * n1**2 + C66 * n3**2, 0, (C13 + C44) * n1 * n3],
        [0, C66 * n1**2 + C44 * n3**2, 0],
        [(C13 + C44) * n1 * n3, 0, C33 * n3**2 + C44 * n1**2],
    ])

    # Solve the eigenvalue problem for phase velocities squared
    eigenvalues, _ = np.linalg.eig(Gamma)
    velocities = np.sqrt(eigenvalues / rho)

    # Compute slowness (s = 1/v)
    slowness = 1 / velocities

    # reorder slowness 
    sorted_indices = np.argsort(slowness)
    slowness_ordered = slowness[sorted_indices]
    
    return slowness_ordered


def phase_velocity_SH(C, rho, theta):
    C11, C12, C13, C33, C44, C66 = C
    n3 = np.cos(np.radians(theta))
    return np.sqrt( (C66*(1-n3**2) + C44*n3**2)/rho )

