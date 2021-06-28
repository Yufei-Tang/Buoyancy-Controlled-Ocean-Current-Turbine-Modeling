# Buoyancy Controlled Ocean Current Turbine Numerical Modeling
This repository contains a numerical modeling of a buoyancy controlled ocean current turbine (OCT), simulated through Simulink in Matlab. The OCT system is a 7 DOF sysetm, including a variable buoyancy tank, variable pitch rotor, main body, and a 607 m mooring cable attached to the ocean floor at a depth of 325 m, as shown in Figure.
<p align="center">
<img src="https://github.com/IRES-FAU/Buoyancy-Controlled-Ocean-Current-Turbine-Modeling/blob/main/Images/OCT_figure-crop_v2.png" width="400">

The **"Numerical Model"** file consists of the following m files:
* ADCP_Init
* attachment_point
* calc_cable_water_vel
* controller_v2
* drag_forces_20m_IBP_v3
* equations_of_motion_cable
* equations_of_motion5_external_generator
* gravity_buoyancy_3_VarBuoy
* rotor_constants_20mD_turbulence
* rotor_forces_20mD_turbulence_v3
* Run_and_test_rotor_model_turbulence
* Single_3d_element_5
* test_Drag_Forces
* turbine_constants_20mD_turbulence_VarBuoy_v2
* turbulence_constants_v3
* waves_spec_calc2
  
The **"Numerical Model"** file also includes a simulink file **Ocean_Current_Turbine_Model_20m_v12_VarBuoy_CL_a** that provides blocks and functions to simulate the numerical simulation of the OCT.

In order to run the codes, you can open and run the simulink (.slx) file.
  
# Citing this repository
Please cite this repository using:

<pre><code>@misc{2021_OCT_IRES,
  Author = {Intelligent and Resilient Energy Systems (IRES) Research Group},
  Doi = {10.5281/zenodo.5035867},
  Howpublished = {https://github.com/IRES-FAU/Buoyancy-Controlled-Ocean-Current-Turbine-Modeling},
  Month = {June},
  Publisher = {Zenodo},
  Title = {Buoyancy Controlled Ocean Current Turbine Numerical Modeling},
  Url = {https://github.com/IRES-FAU/Buoyancy-Controlled-Ocean-Current-Turbine-Modeling},
  Version = {1.0.0},
  Year = {2021}}
</code></pre>

# Pulications
The following publications out of the IRES-FAU research group used/referred to this repository:
