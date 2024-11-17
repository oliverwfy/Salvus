import numpy as np
import xarray as xr

import salvus.namespace as sn
import salvus.mesh.layered_meshing as lm

from Elastic_Material import Austenite

f_c = 5*1e6 

# 2D box domain parameters (length in m)
x_length = 10 * 1e-3
y_length = 20 * 1e-3
x_range = (0., x_length) 
y_range = (0., y_length) 

matl = Austenite()
params = {
    'RHO': matl.RHO,
    'C11': matl.C11,
    'C12': matl.C12,
    'C13': matl.C13,
    'C33': matl.C33,
    'C44': matl.C44
}

mesh = lm.mesh_from_domain(
    domain=sn.domain.dim2.BoxDomain(x0=x_range[0], x1=x_range[1], y0=y_range[0], y1=y_range[1]),
    model=sn.material.elastic.hexagonal.TensorComponents(**params),
    mesh_resolution=sn.MeshResolution(reference_frequency=f_c),
    )


# # define isotropic material with L and S wave Velocities.
# m0 = lm.material.elastic.Velocity.from_params(rho=1.0, vp=1.0, vs=0.5)

# # m0.ds returns the xarray of parameters of m0

# # define isotropic material with Lamé Parameters
# m0_lame = lm.material.elastic.isotropic.LameParameters.from_params(
#     lam=0.5, mu=0.25, rho=1.0
# )
# # m0_lame.ds returns the xarray of parameters of m0

# # all materials can be transformed from one parameterization to another 
# # m0 = lm.material.elastic.Velocity.from_material(m0_lame)

# # compute the minimum wavelength required to determine the required mesh size
# # m0.to_wavelength_oracle(n_dim=3) 
# # print(f"{m0.to_wavelength_oracle(n_dim=3)}\n{m0_lame.to_wavelength_oracle(n_dim=3)}")


# lm_0 = lm.LayeredModel(m0)

# # query the models it contains 
# # lm_0.models 

# # query the interfaces it contains
# # lm_0.interfaces
 
 
# #  If not specified explicitly, Salvus will automatically add bounding hyperplanes to the top and bottom of the model
# print("\n\n".join(str(s) for s in lm_0.complete().interfaces))

# # inspect complete layered model, including bounding interfaces and internal materials,
# print("\n\n".join(str(s) for s in lm_0.complete().strata))

# # The absolute locations of the boundary interfaces are not yet determined
# i0, i1 = lm_0.complete().interfaces
# print(i0.da.reference_elevation)
# print(i1.da.reference_elevation)

# # define domain 2D Box
# d_2d = sn.domain.dim2.BoxDomain(x0=0, x1=1, y0=0, y1=1)
# d_2d.plot()

# # the default order is of 1
# mr = sn.MeshResolution(
#     reference_frequency=2.0, elements_per_wavelength=1.5, model_order=2
# )

# # define layered mesh by domain, model of materials, resolution
# mesh = lm.mesh_from_domain(domain=d_2d, 
#                            model=lm_0, 
#                            mesh_resolution=mr)

# # 3D mesh requires 3D domain  
# d_3d = sn.domain.dim3.BoxDomain(x0=0, x1=1, y0=0, y1=1, z0=0, z1=1)
# mesh_3d = lm.mesh_from_domain(domain=d_3d, model=lm_0, mesh_resolution=mr)

# # the above homogeneous and isotropic elastic model can be solved analytically
# # define a inhomogeneous elastic model with linear velocity gradient in depth 
# vp_grad = xr.DataArray(
#     np.linspace(1.0, 2.0, 11), [("y", np.linspace(1.0, 0.0, 11))]
# )
# m1 = lm.material.elastic.Velocity.from_params(
#     rho=1.0, vp=vp_grad, vs=0.5 * vp_grad
# )

# # 2D mesh of model 1
# mesh_m1 = lm.mesh_from_domain(domain=d_2d, model=m1, mesh_resolution=mr)
# # 3D mesh of model 1
# mesh_3d_m1 = lm.mesh_from_domain(domain=d_3d, model=m1, mesh_resolution=mr)

# # define linear velocity gradient in the z direction for 3D.
# # vp_grad_z = vp_grad.rename({"y": "z"})

# # m1_z = lm.material.elastic.Velocity.from_params(
# #     rho=1.0, vp=vp_grad_z, vs=0.5 * vp_grad_z
# # )

# # lm.mesh_from_domain(domain=d_3d, model=m1_z, mesh_resolution=mr)

# # specify vertical coordinate "v" ("y" for 2D and "z" for 3D)
# vp_grad_v = xr.DataArray(
#     np.linspace(1.0, 2.0, 11), [("v", np.linspace(1.0, 0.0, 11))]
# )
# m1_v = lm.material.elastic.Velocity.from_params(
#     rho=1.0, vp=vp_grad_v, vs=0.5 * vp_grad_v
# )

# # adding layers 
# m2_l0 = lm.material.acoustic.Velocity.from_params(rho=0.5, vp=0.5)
# m2_l1 = m1_v

# # define interface as Hyperplane
# m2_i0 = lm.interface.Hyperplane.at(0.5)

# # initialize layered-model with layers of materials and interface between layers. 
# lm_2 = lm.LayeredModel([m2_l0, m2_i0, m2_l1])

# # unstructured mesh (irregular)
# mesh_2_layers = lm.mesh_from_domain(domain=d_2d, model=lm_2, mesh_resolution=mr)



# # define 2 deformed interfaces
# x = np.linspace(0, 1, 51)

# m3_i0 = lm.interface.Curve.from_points(
#     x,
#     np.sin(2 * np.pi * x) * 0.1 - 0.1,
#     reference_elevation=lm.interface.Depth(0.0),
#     axis="x",
# )

# m3_i1 = lm.interface.Curve.from_points(
#     x,
#     np.sin(np.pi * x) * 0.2,
#     reference_elevation=lm.interface.Depth(0.5),
#     axis="x",
# )

# lm_3 = lm.LayeredModel([m3_i0, m2_l0, m3_i1, m2_l1])

# mesh_deformed_int = lm.mesh_from_domain(domain=d_2d, model=lm_3, mesh_resolution=mr)