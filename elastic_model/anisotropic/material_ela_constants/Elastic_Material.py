import numpy as np
 

# class Elastic_Material:
#   def __init__(self, name):
#     self.name = name
#     if self.name = "Austenite":
      
#   def C
  
  
  
class Austenite:
  def __init__(self):
    self.RHO = 8100
    self.C11 = 217.1 * 1e9
    self.C13 = 144.4 * 1e9
    self.C33 = 263.2 * 1e9
    self.C44 = 82.4 * 1e9
    self.C66 = 128.4 * 1e9
    self.C12 = self.C11 - 2 * self.C66
    self.C = self.tensor()  
    self.params = self.parameters()
    
  def density(self):
    return self.RHO

  def tensor(self, theta=0):
    C = np.zeros((6,6))
    C[0,0] = self.C11
    C[0,1] = self.C12
    C[0,2] = self.C13
    C[1,1] = self.C11
    C[1,2] = self.C33
    C[2,2] = self.C33
    C[3,3] = self.C44
    C[4,4] = self.C44
    C[5,5] = (self.C11 - self.C12) / 2
    
    C = (C + C.T) / 2.0
    return C
  
  def parameters(self, theta=0):
    return {
    'RHO': self.RHO,
    'C11': self.C11,
    'C12': self.C12,
    'C13': self.C13,
    'C33': self.C33,
    'C44': self.C44
    }


def elaTensor_rotation(C, theta=0):
  return C
  
    
