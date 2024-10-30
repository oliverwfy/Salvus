# Full wavefield simulation of Linear Elastic model

[isotropic.py](https://github.com/oliverwfy/Salvus/blob/main/elastic_model/isotropic.py) simulates wave propagation in the isotropic elastic model.

The **domain** is defined as a 2D Box with 10mm x 20mm.

![isotropic_2d_domain](image/isotropic_2d_domain.png)

Then generate mesh and config event (the spatial information of sources and receives), sources in Salvus initialize the force $\mathbf{f}$ in the strong form of elastodynamics equation.

> ## Elastodynamics Equation:
> $$\rho\mathbf{\ddot{u}} = \nabla \cdot \mathbf{C} (\nabla \mathbf{u}) + \mathbf{f}$$ 
><br />

The default mesh alone with a **vector point source** at (5mm, 20mm) and **11 receivers** at the bottom line (**y=0**) is shown below:

![isotropic_2d_mesh](image/isotropic_2d_mesh.png)

The spectral-element method works by approximating the dynamic field variables as **Lagrange polynomials** with time-varying coefficients. These polynomials are defined using the **Gauss-Lobatto-Legendre (GLL)** collocation points defined within each element. Salvus use polynomial of order of degree 4 as default, thus there will be 5 GLL points along each edge of a given element, which results in 25 and 125 points each for 2-D and 3-D elements, respectively.

The **source time function (stf)** used is the Richer wavelet with center frequency at 2MHz and only in the direction of **-y**,

![isotropic_2d_waveform_y](image/isotropic_2d_Ricker.png)


The received **displacement fields** in **y** direction:  

![isotropic_2d_waveform_y](image/isotropic_2d_waveforms_component_y.png)

The received **displacement fields** in **x** direction:  

![isotropic_2d_waveform_x](image/isotropic_2d_waveforms_component_x.png)

Here is an animation of the full wavefield simulation (**magnitude of displacement field**):

![isotropic_2d_animation](image/isotropic_free_surface.gif)
