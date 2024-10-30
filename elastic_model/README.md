# Full wavefield simulation of Linear Elastic model

[isotropic.py](https://github.com/oliverwfy/Salvus/blob/main/elastic_model/isotropic.py) simulates wave propagation in the isotropic elastic model.

The **domain** is defined as a 2D Box with 10mm x 20mm.

![isotropic_2d_domain](image/isotropic_2d_domain.png)

Then config a **vector point source** at (5mm, 20mm) and 11 **receivers** at the bottom line (**y=0**):

![isotropic_2d_mesh](image/isotropic_2d_mesh.png)

The **source time function (stf)** used is the Richer wavelet with center frequency at 2MHz and only in the direction of **-y**:

![isotropic_2d_waveform_y](image/isotropic_2d_Ricker.png)


The received displacement fields in **y** direction:  

![isotropic_2d_waveform_y](image/isotropic_2d_waveforms_component_y.png)

The received displacement fields in **x** direction:  

![isotropic_2d_waveform_x](image/isotropic_2d_waveforms_component_x.png)

Here is an animation for the full wavefield simulation (**magnitude of displacement field**):

![isotropic_2d_animation](image/isotropic_free_surface.gif)
