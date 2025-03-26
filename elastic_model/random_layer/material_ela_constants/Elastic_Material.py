import numpy as np
 




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
    return np.rint(C_rotated)


  
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

    # Calculate epsilon
    epsilon = (C11 - C33) / (2 * C33)

    # Calculate delta
    delta = ((C13 + C44) ** 2 - (C33 - C44) ** 2) / (2 * C33 * (C33 - C44))

    # Calculate eta
    eta = (epsilon - delta) / (1 + 2 * delta)


    # Return results as a dictionary
    return {
        "RHO": rho,
        "VPV": VPV,
        "VPH": VPH,
        "VSV": VSV,
        "VSH": VSH,
        "ETA": eta
    }
    





class Austenite:
    def __init__(self):
        # Define parameters
        self.RHO = 8100  # Density in kg/m^3
        self.C11 = 217.1e9  # Stiffness in Pa
        self.C13 = 144.4e9  # Stiffness in Pa
        self.C33 = 263.2e9  # Stiffness in Pa
        self.C44 = 82.4e9   # Stiffness in Pa
        self.C66 = 128.4e9  # Stiffness in Pa

        # Derived parameters
        self.C12 = self.C11 - 2 * self.C66  # Derived from the relationship
        self.C22 = self.C11                # Symmetry in TTI media
        self.C23 = self.C13                # Symmetry in TTI media
        self.C55 = self.C44                # Symmetry in TTI media
        
        self.C = self.tensor()  

    

    def density(self):
        return self.RHO

    def tensor(self):
        C = np.zeros((6,6))
        C[0,0] = self.C11
        C[0,1] = self.C12
        C[0,2] = self.C13
        C[1,1] = self.C22
        C[1,2] = self.C23
        C[2,2] = self.C33
        C[3,3] = self.C44
        C[4,4] = self.C55
        C[5,5] = (self.C11 - self.C12) / 2
        
        C = (C + C.T) / 2.0
        return C
    
    def rotated_tensor(self, theta):
        return rotated_elasticityTensor(self.tensor(), theta)
        
    def rotated_VTI_approx(self, theta):
        rotated_C = self.rotated_tensor(theta)
        params = {
            "rho": self.RHO, 
            "c11": rotated_C[0, 0],
            "c12": rotated_C[0, 0] - 2 * rotated_C[5, 5],
            "c13": rotated_C[0, 2],
            "c14": 0,
            "c15": 0,
            "c16": 0,
            "c22": rotated_C[0, 0],
            "c23": rotated_C[0, 2],
            "c24": 0,
            "c25": 0,
            "c26": 0,
            "c33": rotated_C[2, 2],
            "c34": 0,
            "c35": 0,
            "c36": 0,
            "c44": rotated_C[3, 3],
            "c45": 0,
            "c46": 0.0,
            "c55": rotated_C[3, 3],
            "c56": 0,
            "c66": rotated_C[5, 5],
        }
        return params

    def rotated_parameters(self, theta):
        rotated_C = self.rotated_tensor(theta)
        return {
        "rho": self.RHO, 
        "c11": rotated_C[0, 0],
        "c12": rotated_C[0, 1],
        "c13": rotated_C[0, 2],
        "c14": rotated_C[0, 3],
        "c15": rotated_C[0, 4],
        "c16": rotated_C[0, 5],
        "c22": rotated_C[1, 1],
        "c23": rotated_C[1, 2],
        "c24": rotated_C[1, 3],
        "c25": rotated_C[1, 4],
        "c26": rotated_C[1, 5],
        "c33": rotated_C[2, 2],
        "c34": rotated_C[2, 3],
        "c35": rotated_C[2, 4],
        "c36": rotated_C[2, 5],
        "c44": rotated_C[3, 3],
        "c45": rotated_C[3, 4],
        "c46": rotated_C[3, 5],
        "c55": rotated_C[4, 4],
        "c56": rotated_C[4, 5],
        "c66": rotated_C[5, 5],
        }






  
class Titanium:
    def __init__(self):
        # Define parameters
        self.RHO = 4506  # Density in kg/m^3
        self.C11 = 162e9  # Stiffness in Pa
        self.C13 = 69e9  # Stiffness in Pa
        self.C33 = 180e9  # Stiffness in Pa
        self.C44 = 47e9   # Stiffness in Pa
        self.C66 = 35e9  # Stiffness in Pa

        # Derived parameters
        self.C12 = self.C11 - 2 * self.C66  # Derived from the relationship
        self.C22 = self.C11                # Symmetry in TTI media
        self.C23 = self.C13                # Symmetry in TTI media
        self.C55 = self.C44                # Symmetry in TTI media
        self.ETA = self.calculate_eta()
        
        self.C = self.tensor()  
        self.params = self.parameters()


    def calculate_eta(self):

        # Calculate epsilon
        epsilon = (self.C11 - self.C33) / (2 * self.C33)

        # Calculate delta
        delta = ((self.C13 + self.C44) ** 2 - (self.C33 - self.C44) ** 2) / (2 * self.C33 * (self.C33 - self.C44))

        # Calculate eta
        eta = (epsilon - delta) / (1 + 2 * delta)
        return eta
    

    def density(self):
        return self.RHO

    def tensor(self):
        C = np.zeros((6,6))
        C[0,0] = self.C11
        C[0,1] = self.C12
        C[0,2] = self.C13
        C[1,1] = self.C22
        C[1,2] = self.C23
        C[2,2] = self.C33
        C[3,3] = self.C44
        C[4,4] = self.C55
        C[5,5] = (self.C11 - self.C12) / 2
        
        C = (C + C.T) / 2.0
        return C
    
    def rotated_tensor(self, theta):
        return rotated_elasticityTensor(self.tensor(), theta)
        
    def parameters(self, theta=None):
        if not theta:
            return {
            'RHO': self.RHO,
            'VPV': np.sqrt(self.C33/self.RHO),
            'VPH': np.sqrt(self.C11/self.RHO),
            'VSV': np.sqrt(self.C55/self.RHO),
            'VSH': np.sqrt(self.C66/self.RHO),
            'ETA': self.ETA}

        else:
            return TTI_velocity_from_tensor(self.rotated_tensor(theta), self.RHO)
        
    def rotated_parameters(self, theta):
        rotated_C = self.rotated_tensor(theta)
        return {
        "C11": rotated_C[0, 0],
        "C12": rotated_C[0, 1],
        "C13": rotated_C[0, 2],
        "C14": rotated_C[0, 3],
        "C15": rotated_C[0, 4],
        "C16": rotated_C[0, 5],
        "C22": rotated_C[1, 1],
        "C23": rotated_C[1, 2],
        "C24": rotated_C[1, 3],
        "C25": rotated_C[1, 4],
        "C26": rotated_C[1, 5],
        "C33": rotated_C[2, 2],
        "C34": rotated_C[2, 3],
        "C35": rotated_C[2, 4],
        "C36": rotated_C[2, 5],
        "C44": rotated_C[3, 3],
        "C45": rotated_C[3, 4],
        "C46": rotated_C[3, 5],
        "C55": rotated_C[4, 4],
        "C56": rotated_C[4, 5],
        "C66": rotated_C[5, 5],
        "RHO": self.RHO,  # Density remains unchanged
        }


