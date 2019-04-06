# Constraint-based Cloth Simulation - Python in Maya

* Implemented a simulation for constraint-based cloth in Maya’s Python API based on Verlet integration. Structural, shear and bend constraints added.
* Added parametric constraints for spheres and planes. 
* Added functionality for saving all keyframes in Maya and rendering a video of the simulation automatically.
* Modular code to allow for easy switching of integrators. All parameters easily accessible. 

Intended use to be by students/researchers, looking to quickly render simulation videos for showcasing.

How to run the code: 
The file is called “clothSim.py”.
1. Paste entire script in Maya’s script editor. A plane mesh (cloth) and a sphere (collider) will be created.
2. Move the sphere to where ever you want it to be. 
3. Open another tab in the script editor and run this command: cloth_sim.time_step(). This command takes the simulation through 1 timestep. Run this as many times.

There is no interface for the the program right now but to turn off the SetKeyFrame functionality, "self.__set_keyframe()" can be commented out in the "time_step(self)" function.

Also, parameters such as the number of constraint solving iterations and time step length can be adjusted by changing the values of the "num_iter" and "t_step" arguments of the "ParticleSystem" function.
