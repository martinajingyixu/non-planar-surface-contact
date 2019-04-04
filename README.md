# 3D-Surface-Contact
## A pipeline to model the 3D frictional surface contact

## 1. Surface contact generation, including:
- parametric surfaces
- discretization of parametric surface
- FEM simulation with two types of soft finger
    - Objects are discribed with non-uniform rational B-Splines (NURBS) for smooth surface generation
    - ANSYS is used to mesh objects and simulate contacts
- **Code**
    - createParametricSurface/
    - FEMSimulation/
- **Required external libraries**
- Create NURBS-based objects
    - MATLAB NURBS Toolbox: https://de.mathworks.com/matlabcentral/fileexchange/26390-nurbs-toolbox-by-d-m-spink
    - MATLAB IGES Toolbox: https://de.mathworks.com/matlabcentral/fileexchange/13253-iges-toolbox
- FEM simulation
    - ANSYS workbench: https://www.ansys.com/
    

## 2. 6D friction computation for non-planar (3D) surfaces
- for parametric surfaces: Numerical integrated wrenches with symbolic integrands
- for meshed surfaces: Sum of the wrench of all elements
- **Code**
    - frictionComputation/
- **Required external libraries**
-  Visualization: show vector with arrow
    - vectarrow: https://de.mathworks.com/matlabcentral/fileexchange/7470-plot-2d-3d-vector-with-arrow?focused=5062305&tab=function

## 3. Fitting limit surface models
### The fitting can apply to both planar and non-planar surfaces with two models:
- Ellipsoidal model
- convex quartic polynomial model
- **Code**
    - fitLimitSurface/
- **Required external libraries**
- convex optimization
    - Convex optimiation cvx: http://cvxr.com/cvx/
- Visualization: show vector with arrow
    - vectarrow: https://de.mathworks.com/matlabcentral/fileexchange/7470-plot-2d-3d-vector-with-arrow?focused=5062305&tab=function


## 4. Mulit-contact modeling
- Sampling the limit surface models
- Build the grasp wrench space
- Check if a 6D wrench disturbance is inside of the grasp wrench space
- **Code**
    - multiContact/
- **Required external libraries**
- check if a disturbance is inside of the qull (grasp stability assessment)
    - inhull: https://de.mathworks.com/matlabcentral/fileexchange/10226-inhull




### License
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.

