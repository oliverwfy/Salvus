# Salvus simulation structure

p (**class**: [Project](https://mondaic.com/docs/2024.1.2/references/python_api/salvus/project/project#project))
- events [list]
    - Event (**class**: [Event](https://mondaic.com/docs/2024.1.2/references/python_api/salvus/flow/collections/event#event) )
        - event_name: [string]
        - sources: [list]
            - VectorPoint2D (**class**: [VectorPoint2D](https://mondaic.com/docs/2024.1.2/references/python_api/salvus/flow/simple_config/source/cartesian#vectorpoint2d) )
                - x [float]
                - y [float]
                - fx [float]
                - fy [float]
        - receivers: [list]
            - Point2D (**class**: [Point2D](https://mondaic.com/docs/2024.1.2/references/python_api/salvus/flow/simple_config/receiver/cartesian#point2d) )
                - x [float]
                - y [float]
                - station_code [string]
                - fields [list[string]]
- sim_config (**class**: [SimulationConfiguration](https://mondaic.com/docs/2024.1.2/references/python_api/salvus/project/configuration/simulation_configuration#simulationconfiguration))
    - name [string]
    - elements_per_wavelength [float]
    - max_frequency_in_hertz [float]
    - model_configuration (**class**: [ModelConfiguration](https://mondaic.com/docs/2024.1.2/references/python_api/salvus/project/configuration/model#modelconfiguration))
        - backgroud_model (**class**: [IsotropicElastic](https://mondaic.com/docs/2024.1.2/references/python_api/salvus/project/configuration/model/background/homogeneous#isotropicelastic)) 
            - rho [float]
            - vp [float]
            - vs [float]
    - absorbing_boundaries (**class**: [AbsorbingBoundaryParameters](https://mondaic.com/docs/2024.1.2/references/python_api/salvus/mesh/simple_mesh/basic_mesh#absorbingboundaryparameters))
        - reference_velocity [float]
        - number_of_wavelengths [float]
        - reference_frequency [float]
        - free_surface [list[string]]
    - event_configuration (**class**: [EventConfiguration](https://mondaic.com/docs/2024.1.2/references/python_api/salvus/project/configuration/event_configuration#eventconfiguration))
        - wavelet (**class**:[Custom](https://mondaic.com/docs/2024.1.2/references/python_api/salvus/flow/simple_config/stf#custom))
        - waveform_simulation_configuration (**class**:[WaveformSimulationConfiguration](https://mondaic.com/docs/2024.1.2/references/python_api/salvus/project/configuration/waveform_simulation_configuration#waveformsimulationconfiguration))
            - start_time_in_seconds [float]
            - end_time_in_seconds [float]
            - spectral_element_order [int]    # default order is 4
            - receiver_sampling_rate_in_time_steps [float]



ed (**class**:[EventData](https://mondaic.com/docs/2024.1.2/references/python_api/salvus/flow/collections/event_data#eventdata))


mesh (**class**:[UnstructuredMesh](https://mondaic.com/docs/2024.1.2/references/python_api/salvus/mesh/data_structures/unstructured_mesh/unstructured_mesh#unstructuredmesh))

CartesianHomogeneousIsotropicElastic2D (**class**:[CartesianHomogeneousIsotropicElastic2D](https://mondaic.com/docs/2024.1.2/references/python_api/salvus/mesh/simple_mesh/basic_mesh#cartesianhomogeneousisotropicelastic2d))
- vp
- vs
- rho 
- x_max
- y_max
- max_frequency 
- elements_per_wavelength 
- tensor_order 
- ab_params (**class**: [AbsorbingBoundaryParameters](https://mondaic.com/docs/2024.1.2/references/python_api/salvus/mesh/simple_mesh/basic_mesh#absorbingboundaryparameters))