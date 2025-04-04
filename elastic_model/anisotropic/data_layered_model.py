import salvus.namespace as sn
from my_code.utilities import *
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import savemat
import h5py
from material_ela_constants.Elastic_Material import Austenite
from scipy.signal import hilbert
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




matl = Austenite() 
rho = matl.RHO
orientation_angle = np.pi/3

c44_ref = matl.C44
c44_rotated = matl.rotated_parameters(orientation_angle)['c44']

v_ref = np.sqrt(c44_ref/rho)
v_rotated= np.sqrt(c44_rotated/rho)



d = 10e-3

# project_name = 'ref_layer'
# simulation_name = "anisotropic_ref_layer"



project_unoriented = 'ref_model'

simulation_unoriented = "mesh_unoriented"

project_interface = 'layered_model'

simulation_interface = "mesh_interface_60"


print("Opening existing project.")
p = sn.Project(path=Path(PROJECT_DIR_WIN, project_unoriented))



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
ed_unoriented = [p.waveforms.get(data_name=simulation_unoriented, events=e)[0] for e in events_list]


time = time_from_ed(ed_unoriented)
time = time.reshape(len(time), -1)

u_y_unoriented = ed_unoriented[0].get_data_cube(receiver_field='displacement', component='Y')[1].T




n_rxs = 1000

u_y_unoriented_mean = [u_y_unoriented[:,i*n_rxs:(i+1)*n_rxs].mean(axis=1) for i in range(5)]



print("Opening existing project.")
p = sn.Project(path=Path(PROJECT_DIR_WIN, project_interface))


ed_interface = [p.waveforms.get(data_name=simulation_interface, events=e)[0] for e in events_list]


time = time_from_ed(ed_interface)
time = time.reshape(len(time), -1)

u_y_interface = ed_interface[0].get_data_cube(receiver_field='displacement', component='Y')[1].T


u_y_interface_mean = [u_y_interface[:,i*n_rxs:(i+1)*n_rxs].mean(axis=1) for i in range(5)]






for i in range(5):
    plt.figure()

    plt.plot(time*1e6, u_y_unoriented_mean[i])
    plt.plot(time*1e6, u_y_interface_mean[i])
    plt.ylabel(rf'$u_2$')
    plt.xlabel('Time (us)')

    plt.legend(['referece', '2 layers'])
    plt.title(fr'Rxs plane $x_3$ ={np.round(10/4*i, 2)} mm')
    plt.savefig(Path(IMAGE_DIR_WIN, fr'com_layered_model_with_interface_x3_{int(10/4*i)}mm.png'))
    # plt.show()





u_mag = np.abs(u_y_unoriented[:,4*n_rxs + int(n_rxs//2)]).reshape(-1,1)
time_start = -0.2*1e-6
time_end = 0.2*1e-6
start_idx = np.abs(time - time_start).argmin()
end_idx = np.abs(time - time_end).argmin()
max_idxs = u_mag[start_idx:end_idx,:].argmax(axis=0) + start_idx
first_time = time[max_idxs].squeeze()



u_y_unoriented_mean_abs = np.abs(u_y_unoriented_mean)
u_y_interface_mean_abs = np.abs(u_y_interface_mean)


i = 4
plt.figure()
plt.plot(time*1e6, u_y_unoriented_mean[i])
plt.plot(time*1e6, u_y_interface_mean[i])
plt.axvline(x = (first_time)*1e6, color='g', linestyle='--', linewidth=1)
plt.axvline(x = (d/v_ref- first_time)*1e6, color='gray', linestyle='--', linewidth=1)
plt.axvline(x = (d/v_ref + d/v_rotated - first_time)*1e6, color='r', linestyle='--', linewidth=1)
plt.axvline(x = (2*d/v_ref - first_time)*1e6, color='b', linestyle='--', linewidth=1)
plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([rf'referece', rf'2 layers'])
plt.title(fr'Rxs plane $x_3$ ={np.round(10/4*i, 2)} mm')
plt.savefig(Path(IMAGE_DIR_WIN, fr'com_layered_model_with_interface_x3_{int(10/4*i)}mm_TOF.png'))




i = 3
plt.figure()
plt.plot(time*1e6, u_y_unoriented_mean[i])
plt.plot(time*1e6, u_y_interface_mean[i])
plt.axvline(x = (d/4/v_ref-first_time)*1e6, color='g', linestyle='--', linewidth=1)
plt.axvline(x = (d/2/v_ref + d/v_rotated + d/4/v_ref - first_time)*1e6, color='r', linestyle='--', linewidth=1)
plt.axvline(x = (7*d/4/v_ref - first_time)*1e6, color='b', linestyle='--', linewidth=1)
plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([rf'referece', rf'2 layers'])
plt.title(fr'Rxs plane $x_3$ ={np.round(10/4*i, 2)} mm')
plt.savefig(Path(IMAGE_DIR_WIN, fr'com_layered_model_with_interface_x3_{int(10/4*i)}mm_TOF.png'))
    
i = 2
plt.figure()
plt.plot(time*1e6, u_y_unoriented_mean[i])
plt.plot(time*1e6, u_y_interface_mean[i])
plt.axvline(x = (d/2/v_ref  - first_time)*1e6, color='g', linestyle='--', linewidth=1)
plt.axvline(x = (d/2/v_ref + d/v_rotated - first_time)*1e6, color='r', linestyle='--', linewidth=1)
plt.axvline(x = (3*d/2/v_ref - first_time)*1e6, color='b', linestyle='--', linewidth=1)
plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([rf'referece', rf'2 layers'])
plt.title(fr'Rxs plane $x_3$ ={np.round(10/4*i, 2)} mm')
plt.savefig(Path(IMAGE_DIR_WIN, fr'com_layered_model_with_interface_x3_{int(10/4*i)}mm_TOF.png'))

i = 1
plt.figure()
plt.plot(time*1e6, u_y_unoriented_mean[i])
plt.plot(time*1e6, u_y_interface_mean[i])
plt.axvline(x = (d/2/v_ref + d/4/v_rotated  - first_time)*1e6, color='g', linestyle='--', linewidth=1)
plt.axvline(x = (d/2/v_ref + d/4/v_ref  - first_time)*1e6, color='gray', linestyle='--', linewidth=1)

plt.axvline(x = (d/2/v_ref + 3*d/4/v_rotated - first_time)*1e6, color='r', linestyle='--', linewidth=1)


plt.axvline(x = (5*d/4/v_ref - first_time)*1e6, color='b', linestyle='--', linewidth=1)


plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([rf'referece', rf'2 layers'])
plt.title(fr'Rxs plane $x_3$ ={np.round(10/4*i, 2)} mm')
plt.savefig(Path(IMAGE_DIR_WIN, fr'com_layered_model_with_interface_x3_{int(10/4*i)}mm_TOF.png'))




i = 0


plt.figure()
plt.plot(time*1e6, u_y_unoriented_mean[i])
plt.plot(time*1e6, u_y_interface_mean[i])
plt.axvline(x = (d/2/v_ref + d/2/v_rotated  - first_time)*1e6, color='g', linestyle='--', linewidth=1)
plt.axvline(x = (d/2/v_ref + d/2/v_ref  - first_time)*1e6, color='gray', linestyle='--', linewidth=1)



plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([rf'referece', rf'2 layers'])
plt.title(fr'Rxs plane $x_3$ ={np.round(10/4*i, 2)} mm')
plt.savefig(Path(IMAGE_DIR_WIN, fr'com_layered_model_with_interface_x3_{int(10/4*i)}mm_TOF.png'))











i = 4
plt.figure()
plt.plot(time*1e6, u_y_unoriented[:,i*n_rxs + int(1/4*n_rxs):i*n_rxs + int(3/4*n_rxs)].mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,i*n_rxs + int(1/4*n_rxs):i*n_rxs + int(3/4*n_rxs)].mean(axis=1))
plt.axvline(x = (first_time)*1e6, color='g', linestyle='--', linewidth=1)
plt.axvline(x = (d/v_ref- first_time)*1e6, color='gray', linestyle='--', linewidth=1)
plt.axvline(x = (d/v_ref + d/v_rotated - first_time)*1e6, color='r', linestyle='--', linewidth=1)
plt.axvline(x = (2*d/v_ref - first_time)*1e6, color='b', linestyle='--', linewidth=1)
plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([rf'referece', rf'2 layers'])
plt.title(fr'Rxs plane $x_3$ ={np.round(10/4*i, 2)} mm')
plt.savefig(Path(IMAGE_DIR_WIN, fr'com_layered_model_with_interface_x3_{int(10/4*i)}mm_TOF_point_rxs.png'))




i = 3
plt.figure()
plt.plot(time*1e6, u_y_unoriented[:,i*n_rxs + int(1/4*n_rxs):i*n_rxs + int(3/4*n_rxs)].mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,i*n_rxs + int(1/4*n_rxs):i*n_rxs + int(3/4*n_rxs)].mean(axis=1))
plt.axvline(x = (d/4/v_ref-first_time)*1e6, color='g', linestyle='--', linewidth=1)
plt.axvline(x = (d/2/v_ref + d/v_rotated + d/4/v_ref - first_time)*1e6, color='r', linestyle='--', linewidth=1)
plt.axvline(x = (7*d/4/v_ref - first_time)*1e6, color='b', linestyle='--', linewidth=1)
plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([rf'referece', rf'2 layers'])
plt.title(fr'Rxs plane $x_3$ ={np.round(10/4*i, 2)} mm')
plt.savefig(Path(IMAGE_DIR_WIN, fr'com_layered_model_with_interface_x3_{int(10/4*i)}mm_TOF_point_rxs.png'))
    
i = 2
plt.figure()
plt.plot(time*1e6, u_y_unoriented[:,i*n_rxs + int(1/4*n_rxs):i*n_rxs + int(3/4*n_rxs)].mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,i*n_rxs + int(1/4*n_rxs):i*n_rxs + int(3/4*n_rxs)].mean(axis=1))
plt.axvline(x = (d/2/v_ref  - first_time)*1e6, color='g', linestyle='--', linewidth=1)
plt.axvline(x = (d/2/v_ref + d/v_rotated - first_time)*1e6, color='r', linestyle='--', linewidth=1)
plt.axvline(x = (3*d/2/v_ref - first_time)*1e6, color='b', linestyle='--', linewidth=1)
plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([rf'referece', rf'2 layers'])
plt.title(fr'Rxs plane $x_3$ ={np.round(10/4*i, 2)} mm')
plt.savefig(Path(IMAGE_DIR_WIN, fr'com_layered_model_with_interface_x3_{int(10/4*i)}mm_TOF_point_rxs.png'))

i = 1
plt.figure()
plt.plot(time*1e6, u_y_unoriented[:,i*n_rxs + int(1/4*n_rxs):i*n_rxs + int(3/4*n_rxs)].mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,i*n_rxs + int(1/4*n_rxs):i*n_rxs + int(3/4*n_rxs)].mean(axis=1))
plt.axvline(x = (d/2/v_ref + d/4/v_rotated  - first_time)*1e6, color='g', linestyle='--', linewidth=1)
plt.axvline(x = (d/2/v_ref + d/4/v_ref  - first_time)*1e6, color='gray', linestyle='--', linewidth=1)

plt.axvline(x = (d/2/v_ref + 3*d/4/v_rotated - first_time)*1e6, color='r', linestyle='--', linewidth=1)


plt.axvline(x = (5*d/4/v_ref - first_time)*1e6, color='b', linestyle='--', linewidth=1)


plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([rf'referece', rf'2 layers'])
plt.title(fr'Rxs plane $x_3$ ={np.round(10/4*i, 2)} mm')
plt.savefig(Path(IMAGE_DIR_WIN, fr'com_layered_model_with_interface_x3_{int(10/4*i)}mm_TOF_point_rxs.png'))




i = 0
plt.figure()
plt.plot(time*1e6, u_y_unoriented[:,i*n_rxs + int(1/4*n_rxs):i*n_rxs + int(3/4*n_rxs)].mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,i*n_rxs + int(1/4*n_rxs):i*n_rxs + int(3/4*n_rxs)].mean(axis=1))
plt.axvline(x = (d/2/v_ref + d/2/v_rotated  - first_time)*1e6, color='g', linestyle='--', linewidth=1)
plt.axvline(x = (d/2/v_ref + d/2/v_ref  - first_time)*1e6, color='gray', linestyle='--', linewidth=1)



plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([rf'referece', rf'2 layers'])
plt.title(fr'Rxs plane $x_3$ ={np.round(10/4*i, 2)} mm')
plt.savefig(Path(IMAGE_DIR_WIN, fr'com_layered_model_with_interface_x3_{int(10/4*i)}mm_TOF_point_rxs.png'))






# plt.savefig(Path(IMAGE_DIR_WIN, fr'com_layered_model_with_interface_bottom_60.png'))


# plt.figure()

# plt.plot(time*1e6, u_y_unoriented_top)
# plt.plot(time*1e6, u_y_interface_top)
# plt.ylabel(rf'$u_y$')
# plt.xlabel('Time (us)')

# plt.legend(['referece (top)', '2 layers (top)'])
# plt.savefig(Path(IMAGE_DIR_WIN, fr'com_layered_model_with_interface_top_60.png'))



















time_start = 0.7*1e-6
time_end = 1.3*1e-6


# for i in range(8):
#     plt.figure()
#     plt.plot(time[:5000]*1e6, u_y[:5000,int(i*45)])
#     plt.axvline(x=d/v_estimated[i*45]*1e6, color='green', linestyle='--', linewidth=1)
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$u_y$')
#     plt.legend(['time signal', 'ToF (SH wave)'])
#     plt.show()
#     # plt.savefig(Path(IMAGE_DIR_WIN, fr'ani_u_y_circle_rxs_{i*45}.png'))


    

# angles_full = np.linspace(0, 360, 360)  # Full circle
# angles_radians = np.radians(angles_full)  # Convert angles to radians for polar plot

# # Plot qP-wave velocity map in polar coordinates
# plt.figure(figsize=(8, 8))
# ax1 = plt.subplot(111, projection='polar')
# ax1.plot(angles_radians, velocity_qP_estimated, label="qP-wave Velocity", color='blue')
# ax1.set_theta_zero_location('N')
# ax1.set_theta_direction(-1)
# ax1.set_title("qP-wave Velocity Polar Plot")
# ax1.legend(loc="upper right")
# # plt.show()
# plt.savefig('qP_wave_Velocity_Polar_Plot')







vpv = 5700.33570594696   
vph = 5177.110114320775
vsv = 3189.4889098682943
vsh = 3981.438414926426





# 2D box domain parameters (length in m)
x_length = 20 * 1e-3
y_length = 20 * 1e-3
x_range = (0., x_length) 
y_range = (0, y_length) 


srcs_pos = [Vector(x_length/2, y_range[1]/2)]
rxs_pos = [Vector(x_length*0.2, y_range[1]),   Vector(x_length/2, y_range[1]),   Vector(x_length*0.8, y_range[1]),
           Vector(x_length*0.2, y_range[1]/2), Vector(x_length/2, y_range[1]/2), Vector(x_length*0.8, y_range[1]/2),
           Vector(x_length*0.2, y_range[0]), Vector(x_length/2, y_range[0]), Vector(x_length*0.8, y_range[0]),]






# idx = 0


# src_pos = srcs_pos[0]
# rx_pos = rxs_pos[idx]
# d_tx_rx = rx_pos.distance(src_pos)

# t_vpv = d_tx_rx/vpv
# t_vsv = d_tx_rx/vsv
# t_vph = d_tx_rx/vph
# t_vsh = d_tx_rx/vsh


# plt.figure()
# plt.plot(time[:5000]*1e6, div_u[:5000,idx])
# plt.xlabel('Time (us)')
# plt.ylabel(rf'$\nabla \cdot u$')
# plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
# plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
# plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
# plt.show()

# plt.figure()
# plt.plot(time[:5000]*1e6, curl_u[:5000,idx])
# plt.xlabel('Time (us)')
# plt.ylabel(rf'$\nabla \times u$')
# plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
# plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
# plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
# plt.show()


# plt.figure()
# plt.plot(time[:5000]*1e6, u_x[:5000,idx])
# plt.xlabel('Time (us)')
# plt.ylabel(rf'$u_x$(m)')
# plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
# plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
# plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
# plt.show()


# plt.figure()
# plt.plot(time[:5000]*1e6, u_y[:5000,idx])
# plt.xlabel('Time (us)')
# plt.ylabel(rf'$u_y$(m)')
# plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
# plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
# plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
# plt.show()






# for idx in range(9):
#     src_pos = srcs_pos[0]
#     rx_pos = rxs_pos[idx]
#     d_tx_rx = rx_pos.distance(src_pos)

#     t_vpv = d_tx_rx/vpv
#     t_vsv = d_tx_rx/vsv
#     t_vph = d_tx_rx/vph
#     t_vsh = d_tx_rx/vsh


#     plt.figure()
#     plt.plot(time[:5000]*1e6, div_u[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$\nabla \cdot u$')
#     plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'ani_middle_grad_u_{idx}.png'))


#     plt.figure()
#     plt.plot(time[:5000]*1e6, curl_u[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$\nabla \times u$')
#     plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'ani_middle_curl_u_{idx}.png'))




#     plt.figure()
#     plt.plot(time[:5000]*1e6, u_x[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$u_x$(m)')
#     plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'ani_middle_displacement_x_{idx}.png'))


#     plt.figure()
#     plt.plot(time[:5000]*1e6, u_y[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$u_y$(m)')
#     plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'ani_middle_displacement_y_{idx}.png'))




# vp = 5700.33570594696
# vs = 3189.4889098682943 
# d = 5*1e-3 

# for idx in range(16):
#     t_vp = d/vp
#     t_vs = d/vs


#     plt.figure()
#     plt.plot(time[:5000]*1e6, div_u[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$\nabla \cdot u$')
#     plt.axvline(x=t_vp*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vs*1e6, color='red', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VP)', 'ToF(VS)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'iso_middle_divergence_u_{idx}.png'))


#     plt.figure()
#     plt.plot(time[:5000]*1e6, curl_u[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$\nabla \times u$')
#     plt.axvline(x=t_vp*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vs*1e6, color='red', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VP)', 'ToF(VS)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'iso_middle_curl_u_{idx}.png'))




#     plt.figure()
#     plt.plot(time[:5000]*1e6, u_x[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$u_x$(m)')
#     plt.axvline(x=t_vp*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vs*1e6, color='red', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VP)', 'ToF(VS)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'iso_middle_displacement_x_{idx}.png'))


#     plt.figure()
#     plt.plot(time[:5000]*1e6, u_y[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$u_y$(m)')
#     plt.axvline(x=t_vp*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vs*1e6, color='red', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VP)', 'ToF(VS)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'iso_middle_displacement_y_{idx}.png'))








# fmc_data, time, rxs_loc, srcs_loc = fmc_data_from_ed(event_data=ed, save_dir=DATA_DIR)

# # save fmc data as .mat
# fmc = {
#     'sim_data': {
#         'time_data': fmc_data,
#         'time':time,
#         'rxs_loc':rxs_loc,
#         'srcs_loc':srcs_loc   
#     }
# }

# mfile_name = 'fmc.mat'

# # savemat(Path(DATA_DIR_WIN, mfile_name), fmc)

# # get all recerived data from #m transducer
# id_m = 1

# # plt.imshow(fmc_data[:, id_m:(id_m+1)*len(rxs_loc)].T, aspect='auto', extent=[time[0]*1e6, time[-1]*1e6, 0, len(rxs_loc)], cmap='viridis')
# # plt.colorbar(label=r'$u_y$')
# # plt.xlabel('Time (us)')
# # plt.ylabel('# receiver')
# # plt.title(f'FMC Data in #{0} src')
# # plt.savefig(Path(IMAGE_DIR, 'fmc_first_src.png'))





