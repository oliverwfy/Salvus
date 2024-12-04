import matplotlib.pyplot as plt
from material_ela_constants.Elastic_Material import Austenite
import numpy as np
from pathlib import Path


IMAGE_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/image'


matl = Austenite()      # load material's elasticity tensor


VPV = matl.params['VPV']
VPH = matl.params['VPH']
VSV = matl.params['VSV']
VSH = matl.params['VSH']

angles = np.linspace(0, 360, 360)
VP_values = []
VS_values = []

rads = angles/180*np.pi


for theta in rads:
    paras = matl.parameters(theta=theta)
    VP_values.append(paras['VPV'])
    VS_values.append(paras['VSV'])





# Polar plot for P-wave velocities
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6, 6))

# Plot P-wave velocity
ax.plot(rads, VP_values, label='P-wave Velocity', color='blue')
ax.set_title('P-wave Velocity (TTI)', va='bottom')
ax.set_theta_zero_location('N')  # 0 degrees at North
ax.set_theta_direction(-1)  # Clockwise

# Annotate VPV and VPH
ax.annotate(f'VPV ({VPV:.0f})', xy=(0, VPV), xytext=(0.1, VPV + 200),
            arrowprops=dict(facecolor='black', arrowstyle='->'))
ax.annotate(f'VPH ({VPH:.0f})', xy=(np.pi / 2, VPH), xytext=(np.pi / 2 + 0.1, VPH + 200),
            arrowprops=dict(facecolor='black', arrowstyle='->'))

# Save and show plot
plt.tight_layout()
plt.savefig(Path(IMAGE_DIR_WIN, 'p_wave_velocity_map.png'))  # Save to specified path
plt.show()

# Show plots
plt.tight_layout()
plt.savefig(Path(IMAGE_DIR_WIN, fr'velocity_map.png'))
plt.show()
