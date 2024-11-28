import numpy as np
 

  
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
        self.ETA = self.calculate_eta()
        
        self.C = self.tensor()  
        self.params = self.parameters()


    def calculate_eta(self):

      eta = ((self.C11 - self.C33) / (2 * self.C33)) - \
            ((self.C13 + self.C44)**2 - (self.C33 - self.C44)**2) / \
            (2 * self.C33 * (self.C33 - self.C44))
      return eta
    

    def density(self):
        return self.RHO

    def tensor(self, theta=0):
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
  
    def parameters(self, theta=0):
      
        return {
        'RHO': self.RHO,
        'VPV': np.sqrt(self.C33/self.RHO),
        'VPH': np.sqrt(self.C11/self.RHO),
        'VSV': np.sqrt(self.C55/self.RHO),
        'VSH': np.sqrt(self.C66/self.RHO),
        'ETA': self.ETA}

  

