import matplotlib.pyplot as plt
from material_ela_constants.Elastic_Material import *
import numpy as np
from my_code.utilities import *
from pathlib import Path


IMAGE_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/image'


matl = Austenite()      # load material's elasticity tensor

C_matrix = matl.C
C11, C12, C13, C33, C44, C66 = C_matrix[0,0], C_matrix[0,1], C_matrix[0,2],\
    C_matrix[2,2], C_matrix[3,3], C_matrix[5,5]

C = [C11, C12, C13, C33, C44, C66]


rho = matl.RHO  



# # Plot qP-wave slowness
# plt.figure(figsize=(8, 8))
# ax1 = plt.subplot(111, projection='polar')
# ax1.plot(angles_radians, slowness_qP, label="qP-wave Slowness", color='green')
# ax1.set_theta_zero_location('N')
# ax1.set_theta_direction(-1)
# ax1.set_title("qP-wave Slowness Polar Plot")
# ax1.legend(loc="upper right")
# plt.show()

# # Plot qSV-wave slowness
# plt.figure(figsize=(8, 8))
# ax2 = plt.subplot(111, projection='polar')
# ax2.plot(angles_radians, slowness_qSV, label="qSV-wave Slowness", color='red')
# ax2.set_theta_zero_location('N')
# ax2.set_theta_direction(-1)
# ax2.set_title("qSV-wave Slowness Polar Plot")
# ax2.legend(loc="upper right")
# plt.show()

# # Plot SH-wave slowness
# plt.figure(figsize=(8, 8))
# ax3 = plt.subplot(111, projection='polar')
# ax3.plot(angles_radians, slowness_SH, label="SH-wave Slowness", color='blue')
# ax3.set_theta_zero_location('N')
# ax3.set_theta_direction(-1)
# ax3.set_title("SH-wave Slowness Polar Plot")
# ax3.legend(loc="upper right")
# plt.show()


# Compute velocity from slowness (v = 1 / s)
# velocity_qSV = 1 / slowness_qSV
# velocity_SH = 1 / slowness_SH

# # Plot qP-wave velocity map in polar coordinates
# plt.figure(figsize=(8, 8))
# ax1 = plt.subplot(111, projection='polar')
# ax1.plot(angles_radians, velocity_qP, label="qP-wave Velocity", color='blue')
# ax1.set_theta_zero_location('N')
# ax1.set_theta_direction(-1)
# ax1.set_title("qP-wave Velocity Polar Plot")
# ax1.legend(loc="upper right")
# # plt.show()
# plt.savefig('qP_wave_Velocity_Polar_Plot')

# # # Plot qSV-wave velocity map in polar coordinates
# # plt.figure(figsize=(8, 8))
# # ax2 = plt.subplot(111, projection='polar')
# # ax2.plot(angles_radians, velocity_qSV, label="qSV-wave Velocity", color='red')
# # ax2.set_theta_zero_location('N')
# # ax2.set_theta_direction(-1)
# # ax2.set_title("qSV-wave Velocity Polar Plot")
# # ax2.legend(loc="upper right")
# # plt.show()

# # # Plot SH-wave velocity map in polar coordinates
# # plt.figure(figsize=(8, 8))
# # ax3 = plt.subplot(111, projection='polar')
# # ax3.plot(angles_radians, velocity_SH, label="SH-wave Velocity", color='green')
# # ax3.set_theta_zero_location('N')
# # ax3.set_theta_direction(-1)
# # ax3.set_title("SH-wave Velocity Polar Plot")
# # ax3.legend(loc="upper right")
# # plt.show()












import salvus.namespace as sn
from my_code.utilities import *
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import savemat
import h5py

# Directories in WSL
PROJECT_DIR = '/home/oliver/workspace/Salvus/elastic_model/anisotropic/Project'
IMAGE_DIR = '/home/oliver/workspace/Salvus/elastic_model/anisotropic/image'
DATA_DIR = '/home/oliver/workspace/Salvus/elastic_model/anisotropic/data'

# Directories in Windows
PROJECT_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/Project'
DATA_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/data'
IMAGE_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/image'

# create dir if it does not exist
Path(IMAGE_DIR).mkdir(parents=True, exist_ok=True)
Path(DATA_DIR).mkdir(parents=True, exist_ok=True)
Path(IMAGE_DIR_WIN).mkdir(parents=True, exist_ok=True)
Path(DATA_DIR_WIN).mkdir(parents=True, exist_ok=True)




# project_name = 'ref_layer'
# simulation_name = "anisotropic_ref_layer"


project_name = 'rotated_triclinic'
simulation_name = "mesh_unoriented"


print("Opening existing project.")
p = sn.Project(path=Path(PROJECT_DIR_WIN, project_name))



# p.viz.nb.waveforms("fmc_simulation", receiver_field="displacement")

# # displacement in y direction
# p.waveforms.get(data_name="fmc_simulation", events=["event_0"])[0].plot(
#     component="Y", receiver_field="displacement"
# )


# # displacement in x direction
# p.waveforms.get(data_name="fmc_simulation", events=["event_0"])[0].plot(
#     component="X", receiver_field="displacement"
# )


# get events from project in correct order 
events_list = reorder_events_list(p.events.list())
ed = [p.waveforms.get(data_name=simulation_name, events=e)[0] for e in events_list]


time = time_from_ed(ed)
time = time.reshape(len(time), -1)



# grad_u = {
#     "0": ed[0].get_data_cube(receiver_field='gradient-of-displacement', component='0')[1].T,
#     "1": ed[0].get_data_cube(receiver_field='gradient-of-displacement', component='1')[1].T,
#     "2": ed[0].get_data_cube(receiver_field='gradient-of-displacement', component='2')[1].T,
#     "3": ed[0].get_data_cube(receiver_field='gradient-of-displacement', component='3')[1].T,
# }

# div_u = grad_u['0'] + grad_u['3']
# curl_u = grad_u['2'] - grad_u['1']


u_x_mid = ed[1].get_data_cube(receiver_field='displacement', component='X')[1].T
u_y_mid = ed[1].get_data_cube(receiver_field='displacement', component='Y')[1].T
u_mag_mid = np.sqrt(u_x_mid**2 + u_y_mid**2)

time_start = -0.3*1e-6
time_end = 0.3*1e-6
start_idx = np.abs(time - time_start).argmin()
end_idx = np.abs(time - time_end).argmin()
max_idxs = u_mag_mid[start_idx:end_idx,:].argmax(axis=0) + start_idx

first_time = time[max_idxs]

u_x = ed[0].get_data_cube(receiver_field='displacement', component='X')[1].T
u_y = ed[0].get_data_cube(receiver_field='displacement', component='Y')[1].T
u_y_abs = np.abs(u_y)
u_mag = np.sqrt(u_x**2 + u_y**2)
# find the FOT in the first arriving envolop

time_start = 0.7*1e-6
time_end = 1.3*1e-6
d = 4*1e-3 

start_idx = np.abs(time - time_start).argmin()
end_idx = np.abs(time - time_end).argmin()



max_idxs = u_y_abs[start_idx:end_idx,:].argmax(axis=0) + start_idx

v_estimated = d/(time[max_idxs]-first_time)



new_order = list(range(270, -1, -1))
new_order += list(range(359, 270, -1))

vp = 5700.33570594696


# # Plot qP-wave velocity map in polar coordinates
# plt.figure(figsize=(8, 8))
# ax1 = plt.subplot(111, projection='polar')
# ax1.plot(angles_radians, [vp]*len(angles_radians), label="P-wave Velocity", color='blue')
# ax1.plot(angles_radians, v_estimated[new_order], label="Estimated Group Velocity", color='red')

# ax1.set_theta_zero_location('N')
# ax1.set_theta_direction(-1)
# ax1.set_title("P-wave Velocity Polar Plot")
# ax1.legend(loc="upper right")
# plt.savefig(Path(IMAGE_DIR_WIN, 'comparison_wave_Velocity_Polar_Plot_x_dir_iso.png'))



# Convert slowness data to polar coordinates
angles_full = np.linspace(0, 360, 360)  # Full circle
slowness_results_full = np.array([compute_slowness(C, rho, theta) for theta in angles_full])

# Extract slowness values for the three wave modes
slowness_qP = slowness_results_full[:, 0]  # quasi-P wave
slowness_qSV = slowness_results_full[:, 1]  # quasi-SV wave
slowness_SH = slowness_results_full[:, 2]  # SH wave


angles_radians = np.radians(angles_full)  # Convert angles to radians for polar plot

v_sh = np.array([phase_velocity_SH(C,rho, theta) for theta in angles_full])

rotate_angle = 90
v_sh_rotated = np.concatenate( (v_sh[-rotate_angle:],v_sh[:-rotate_angle]) ) 

rotated_para = matl.rotated_parameters(np.radians(rotate_angle))

a = rotated_para['c44']
b = rotated_para['c66']

v_rot = np.sqrt( (a*np.cos(angles_radians)**2 + b*np.sin(angles_radians)**2)/rho )


# Plot qP-wave velocity map in polar coordinates
plt.figure(figsize=(8, 8))
ax1 = plt.subplot(111, projection='polar')
ax1.plot(angles_radians, v_sh_rotated, label="SH Velocity map", color='blue')
ax1.set_theta_zero_location('N')
ax1.set_theta_direction(-1)
ax1.plot(angles_radians, v_rot, 'r--', label="Approximation")
ax1.set_title("SH wave Velocity Polar Plot")
ax1.legend(loc="upper right")
plt.savefig(Path(IMAGE_DIR_WIN, f'approxiamtion_SH_Velocity_Polar_Plot_y_dir_rotate_{rotate_angle}.png'))



# # Plot qP-wave velocity map in polar coordinates
# plt.figure(figsize=(8, 8))
# ax1 = plt.subplot(111, projection='polar')
# ax1.plot(angles_radians, v_sh_rotated, label="SH Velocity map", color='blue')
# ax1.set_theta_zero_location('N')
# ax1.set_theta_direction(-1)
# ax1.plot(angles_radians, v_estimated, label="Estimated Group Velocity", color='red')
# ax1.set_title("SH wave Velocity Polar Plot")
# ax1.legend(loc="upper right")
# # plt.show()
# plt.savefig(Path(IMAGE_DIR_WIN, f'approxiamtion_SH_Velocity_Polar_Plot_y_dir_rotate_{rotate_angle}.png'))



# plt.savefig(Path(IMAGE_DIR_WIN, 'comparison_SH_Velocity_Polar_Plot_y_dir.png'))

# # Plot qP-wave velocity map in polar coordinates
# plt.figure(figsize=(8, 8))
# ax1 = plt.subplot(111, projection='polar')
# ax1.plot(angles_radians, 1/np.array(v_sh), label="SH Velocity map", color='blue')
# ax1.set_theta_zero_location('N')
# ax1.set_theta_direction(-1)
# ax1.set_title("SH wave slowness Polar Plot")
# plt.show()