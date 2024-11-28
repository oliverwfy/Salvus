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
    fmc_data = np.zeros([len(time), len(srcs_loc)*len(rxs_loc)])

    for i, _ in enumerate(event_data):
        fmc_data[:,i*len(rxs_loc):(i+1)*len(rxs_loc)] = \
        event_data[i].get_waveform_data_xarray('displacement').sel(components="Y").values.T
    
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
        




  
def rotation_matrix_y(theta):
    """
    Computes the 3D rotation matrix about the y-axis for a given angle theta.

    Parameters:
    theta (float): Rotation angle in radians.

    Returns:
    numpy.ndarray: 3x3 rotation matrix.
    """
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    return np.array([
        [cos_theta,  0, sin_theta],
        [0,          1, 0],
        [-sin_theta, 0, cos_theta]
    ])


def transformation_T(Q):
    """
    Computes the transformation matrix T (6x6) in Voigt notation
    from a 3D rotation matrix Q.

    Parameters:
    Q (numpy.ndarray): 3x3 rotation matrix.

    Returns:
    numpy.ndarray: 6x6 transformation matrix.
    """
    T = np.array([
        [Q[0, 0]**2, Q[0, 1]**2, Q[0, 2]**2, 2*Q[0, 1]*Q[0, 2], 2*Q[0, 0]*Q[0, 2], 2*Q[0, 0]*Q[0, 1]],
        [Q[1, 0]**2, Q[1, 1]**2, Q[1, 2]**2, 2*Q[1, 1]*Q[1, 2], 2*Q[1, 0]*Q[1, 2], 2*Q[1, 0]*Q[1, 1]],
        [Q[2, 0]**2, Q[2, 1]**2, Q[2, 2]**2, 2*Q[2, 1]*Q[2, 2], 2*Q[2, 0]*Q[2, 2], 2*Q[2, 0]*Q[2, 1]],
        [2*Q[1, 0]*Q[2, 0], 2*Q[1, 1]*Q[2, 1], 2*Q[1, 2]*Q[2, 2], Q[1, 1]*Q[2, 2] + Q[1, 2]*Q[2, 1], Q[1, 0]*Q[2, 2] + Q[1, 2]*Q[2, 0], Q[1, 0]*Q[2, 1] + Q[1, 1]*Q[2, 0]],
        [2*Q[0, 0]*Q[2, 0], 2*Q[0, 1]*Q[2, 1], 2*Q[0, 2]*Q[2, 2], Q[0, 1]*Q[2, 2] + Q[0, 2]*Q[2, 1], Q[0, 0]*Q[2, 2] + Q[0, 2]*Q[2, 0], Q[0, 0]*Q[2, 1] + Q[0, 1]*Q[2, 0]],
        [2*Q[0, 0]*Q[1, 0], 2*Q[0, 1]*Q[1, 1], 2*Q[0, 2]*Q[1, 2], Q[0, 1]*Q[1, 2] + Q[0, 2]*Q[1, 1], Q[0, 0]*Q[1, 2] + Q[0, 2]*Q[1, 0], Q[0, 0]*Q[1, 1] + Q[0, 1]*Q[1, 0]]
    ])
    return T


def rotated_elasticityTensor(C, theta):
    """
    Rotates the elasticity tensor C (6x6 in Voigt notation) by angle theta
    about the y-axis.

    Parameters:
    C (numpy.ndarray): 6x6 elasticity tensor in Voigt notation.
    theta (float): Rotation angle in radians.

    Returns:
    numpy.ndarray: Rotated 6x6 elasticity tensor in Voigt notation.
    """
    Q = rotation_matrix_y(theta)  # Compute the rotation matrix
    T = transformation_T(Q)      # Compute the transformation matrix (6x6)
    C_rotated = T @ C @ T.T       # Perform the tensor rotation in Voigt notation
    return C_rotated
  
  
def TTI_velocity_from_tensor(C, rho):
    """
    Computes velocities (VPV, VPH, VSV, VSH) and ETA (η) from a 3D TTI stiffness tensor.

    Parameters:
    C (numpy.ndarray): 6x6 stiffness tensor in Voigt notation.
    rho (float): Density of the material (kg/m^3).

    Returns:
    dict: Dictionary containing VPV, VPH, VSV, and ETA.
    """
    # Extract relevant stiffness coefficients
    C11 = C[0, 0]  # Horizontal P-wave stiffness
    C33 = C[2, 2]  # Vertical P-wave stiffness
    C44 = C[3, 3]  # Vertical S-wave stiffness
    C55 = C44      # Symmetry assumption: C55 = C44
    C66 = C[5,5]   # Horizontal S-wave stiffness
    C13 = C[0, 2]  # Coupling stiffness

    # Calculate velocities
    VPV = np.sqrt(C33 / rho)  # Vertical P-wave velocity
    VPH = np.sqrt(C11 / rho)  # Horizontal P-wave velocity
    VSV = np.sqrt(C55 / rho)  # Vertical S-wave velocity
    VSH = np.sqrt(C66 / rho)  # Horizontal S-wave velocity

    # Calculate ETA
    eta_numerator = (C11 - C33) / (2 * C33)
    eta_denominator = ((C13 + C44)**2 - (C33 - C44)**2) / (2 * C33 * (C33 - C44))
    eta = eta_numerator - eta_denominator

    # Return results as a dictionary
    return {
        "RHO": rho,
        "VPV": VPV,
        "VPH": VPH,
        "VSV": VSV,
        "VSH": VSH,
        "ETA": eta
    }    