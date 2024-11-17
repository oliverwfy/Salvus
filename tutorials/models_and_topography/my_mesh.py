import salvus.namespace as sn




class Isotropic_Material:
    def __init__(self, RHO, VP, VS):
        self.VS = VS
        self.VP = VP
        self.RHO = RHO


def get_basic_mesh(matlt: Isotropic_Material, x_max: float, y_max: float
                   , centre_freq: float, elements_per_wavelength = 2, 
                   tensor_order=4, abs_bdry=None) -> sn.UnstructuredMesh:
    
    max_frequency = 2 * centre_freq        

    mesh = sn.simple_mesh.basic_mesh.CartesianHomogeneousIsotropicElastic2D(
        vp=matlt.VP, vs=matlt.VS, rho=matlt.RHO, x_max=x_max, y_max = y_max, 
        max_frequency = max_frequency, elements_per_wavelength =elements_per_wavelength,  
        tensor_order = tensor_order, 
        ab_params = abs_bdry    
    ).create_mesh()

    return mesh     