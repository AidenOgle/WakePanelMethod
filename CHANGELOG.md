# CHANGELOG

All notable changes to this project will be documented in this file

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - YYYY/MM/DD

This update primarily addresses optimization of vortex influence and wake roll-up through vectorization of nested loops, reducing time complexity from $`O(nm^2+n^2m+n^3)`$ to $`O(n)`$. Although this new time complexity has not been rigourously verified, the difference in program run time is noticably significant and allows for simulation lengths that were previously infeasable.

In the process of restructuring for vectorization, sections were reformatted to provide clarity on underlying mechanisms. Notably, the hardcoded vortex core size correction in the Wake Roll-Up section was changed to a function to allow for user-definable core sizes and Vortex Blob Regularization methods. This is accompanied by functions that model vortex spreading and dissipation respectively to simulate aspects of visocity. Additionally, setup functions were modularized to `initialization.py` and their input parameters were concentrated to improve readability and the ease of user-input.

### Added

- Implemented Vortex Blob Method to smooth vortex singularities and allow for core-size dependent visocus modeling. Options included for both Cauchy `'linear'` and Lamb-Oseen `'gaussian'` methods. Replaced built-in linear blob in wake roll-up section
    - Blob regularization is only applied to wake vortices to avoid boundary condition degradation caused by applying core-sizes to bound body vortices. Initial core-size is adjusted by the `blob_kernel` variable scaled against a user defined length scale
    
- Implemented Wake Vortex Spreading mechanic option to simulate viscous effects. Vortex core size $`\delta(t)`$ increases with $`\sqrt{4 \nu t}`$ (kinematic viscosity dependency). Requires a Vortex Blob Method to be active to have an effect
  
- Implemented Wake Vortex Decay mechanic option to simulate viscous effects. Options included for both Exponential Timescale decay and Diffusive decay. Diffusive Decay method requires vortex spreading to be active via `spreading = True`
    - `'exponential'` is scaled by a time constant `tau` dimensionalized by `simlength`, whereas `'diffusive'` is scaled by an effective radius `r_eff` dimensionalized by average panel length. Both are scalable by the term `decay_kernel`
      
- Added `drag` and `coeff_drag` terms as a result of decomposition the lift terms into the global reference frame. Both `lift` and `drag` calculations are experimental in this version because it is unclear if the underlying theories are being applied correctly. The calculation of forces on the body and their related coefficients needs significant work to ensure their implementation is consistent with theory and the assumptions being made

- Created module `initialization.py` to contain `generateNACAfoil(...)`, `body_creation(...)`, `pitching(...)`, and `vortexblobinit(...)`. In the future could add flow velocity initialization as a function with a similar structure to `pitching(...)`

### Changed

- Fully vectorized the Induced Velocity Coefficient matrices to avoid repeatedly iterating through panels and wakepoints
  
- Fully vectorized the Wake Roll-Up section to avoid doubly iterating over all wakepoints
    - Eliminated self-influence terms from wake-wake interactions in the Wake Roll-Up section
      
- Changed the plotting of vortex points to represent their local circulation intensity using the same color scale as the wake vortices
  
- When `spreading = True`, wakepoints are plotted for their marker size to grow proportionally to their percent change in core area. Marker size is not directly physically meaningful and is scaled by the factor `spreading_displayscale` to highlight the normally imperceptable change in vortex core size. To disable this effect, set `spreading_displayscale` to zero
  
- Added `'bodyshape'` variable to be assigned during the body creation to flag bodies as either `'open'` or `'closed'`. This flag is then used in the enforcement of the Kutta Condition
  
- Kutta Condition enforcement in the creation of A and B matricies was changed to support the geometry of closed bodies. The open body (original) method enforced a global circulation constraint such that total circulation is conserved between the body and wake. The new closed body method directly enforces the strong version of the Kutta condition by prescribing the upper and lower trailing edge vorticies to have the same finite circulation
    - Current code formulation applies the Kutta condition regardless of body shape. This results in physically unrealistic flows for bluff bodies like the `'circle'`. Intend to make the enforcement of the Kutta condition toggleable a future update and allow dynaimic definition of separation points
      
- Changed the `lift` and `coeff_lift` calculations to account for pitching and angular components. Decomposed the steady and unsteady terms into their lift and drag components in the global reference frame

### Fixed

- Fixed aspect ratio distortion of the plot resulting from the addition of the color bar
  
- Changed plotting to GridSpec layout to fix plot formatting issues and allow for easier adjustment. Added plot sizing parameters
  
- Fixed the `unsteady` term to account for nonuniform panel sizing

## [0.2.0] - 2025/06/15

### Added

- FLUID SIMULATION: Incorporated vortex panel method used in original MATLAB code to find lift derived from body circulation and wake behavior. The main difference between the implementations in MATLAB and here is that the MATLAB implementation has local variables that are plotted and written over for each timestep iteration, whereas the implementation here stores values outside of the loop for each time step. This change leads to slightly decreased performance but better plotting stability and access to data at different time steps
  
  - BODY SIMULATION:
    - Have to double check if Lift (L) and the Coefficient of Lift (CL) are calculated correctly in accordance with nondimensionalization, as the original MATLAB version didn't allow for variable velocity or density
      
    - Intend to create animated plots for L and CL and generate data for other properties like drag
      
  - WAKE SIMULATION:
      - the MATLAB implementation structured the assignment of values for wake points as a queue, with the newest point generated by the trailing edge being the first index. I believe this is done so that in the Wake Roll-Up section, flow information is propagated backwards along the wake during iteration to be consistent with the physical process. This methodology is applied in this version as well, but this requires an obtuse indexing system in the plotting section
        
      - Intend to implement prescribed viscocity effect for wake circulation dissipation
        
  - COMPLEXITY and OPTIMIZATION:
    - The total complexity of the program is $`O(nm^2+n^2m+n^3)`$, where 'n' is the number of timesteps and 'm' is the number of panels. In most scenarios, timesteps will dominate over the number of panels, resulting in a complexity of $`O(n^3)`$.
      
    - Optimization would ideally be accomplished through vectorization to eliminate extraneous loops. This should be possible for both the Wake Roll-Up and Induced Velocity Coefficient sections, but will require some refamiliarization. I will explore this option in the next update
      
    - Intend to implement parameters that allow wake calculation to be turned off or reduced. Reduction could be in the form of  time resolution by only generating/calculating wake points every nth time step, or by having a range of influence such that wake points are only compared against points that are close in index (i.e. +/-50 index)
  
- NACA AIRFOIL : Added a 4-Digit NACA airfoil option for body geometry to directly generate body profiles based on the 4-digit series airfoils. Supports both symmetric and cambered airfoils. [The equations on this site](http://airfoiltools.com/airfoil/naca4digit) were used as reference. This is an experimental feature, it may be more reliable to directly input x-y points using the `'custom'` method

  - Incorporates new argument `trailingedge=` which allows the user to specifiy if the trailing edge has a finite thickness or is closed. Viable assignments are `'open'`, `'closed'`, and `'prescribed'`. `'prescribed'` and `'closed'` should achieve the same thing, but `'closed'` uses an equation whereas `'prescribed'` directly assigns the thickness at the trailing edge to be zero. This assignment may be depreciated in the future depending on usefulness
    
  - Due to the way vortex and collocation points are defined relative to body points based on order, the top and bottom surfaces have an offset in the order of vortex and collocation points compared to eachother. This is most noticeable for symmetric airfoils with a low number of panels (e.g. `NACA='0015', num_panels=10, display_vortex=True, display_colloc=True`)
    
  - Odd numbers of panels would cause a duplicate point at either the leading or trailing edges, resulting in coincident vortex and collocation points. To avoid this, the number of panels is shunted up to an even value. If an odd number of panels is desired, it is reccomended to use the 'custom' method  to manually define points

- CUSTOM BODY : Option to define custom body geometry using x-y coordinates. Currently inputted as two arrays, may change to have different syntax in the future to reflect how other programs export x-y data

- CUSTOM PITCHING : Option to define custom pitching profile using an array of values corresponding to timesteps. Allows for pitching profiles to be defined by functions/equations beyond the preconstructed ones provided

- Added a progress bar and program timer. Unlike the original MATLAB code, the results are not shown as they are generated and are instead displayed asynchronously after the program has finished. This setup improves the viewability of the results animation but leads to a longer wait time between initialization and output, especially for large panel numbers. The progress bar and program timer were implemented for quality of life and future optimization purposes

### Changed

- Reformatted the geometry creation section into a function structure as `body_creation`
  
- Changed how arguments are passed to the `pitching` function to be more consistent. Now uses a keyword argument structure
  
- Added an exception handling case for improper argument syntax and inputs for the `pitching` and `body_creation` functions to account for the user-unfriendliness of the keyword argument structure. Might implement prompt system in the future for easier user-input
  
- Added `simlength` variable to drive `num_step` to better reflect the time in whole units instead of hundreths. Allows floats. No direct usage yet, added primarily as a reference variable
  
### Fixed

- Fixed `'periodic'` case for `pitching` function to reflect direction conventions. The positive anglular direction for the body reference frame cooresponds to the negative global anglular direction from cartesian frame, which was not originally implemented correctly. Changed `phase` argument to default to `0` to reflect intended default motion under currected system.

## [0.1.0] -  2025/06/04
Created initial version of Panel Method code
Only handles geometry, motion (displacement & pitching), and plotting animation so far. Wake shedding and pressure profile will be addressed in next update

### Added

- GEOMETRY:
  - created presets for `'line'`, `'ellipse'`, and `'circle'` bodies
  - currently bodies are initially centered at the origin (either about their center or quarter-chord). Plan on allowing initial offsets in the future
- DISPLACEMENT:
  - changed from MATLAB code to allow for irregular displacements. Uses stepwise motion
  - currently only allows for only constant velocity or linear accel/deceleration
  - planned to change to follow same structure as PITCHING section and allow for more predefined or custom velocity profiles
- PITCHING:
  - added predefined pitching profiles `'constant'`, `'constant rate between values'`, `'constant rate using rate'`, and `'periodic'`
  - currently structured as a function that passes `args`, will be changed in the future to be more consistent between cases
- PLOTTING:
  - uses FuncAnimation module from matplotlib library
  - added option to automatically scale x-axis to fit the displacement of the body
  - added plotscaling calculation to keep consistent scale between axes and prevent graphical distortion of the body
    - if the auto axis scaling causes the plotsize to exceed the display size, the scaling becomes inconsistent. Will fix next update
  -added options to display/hide the body, vortex points, and collocation points respectively


[0.2.0]: https://github.com/AidenOgle/2DPanelMethod/blob/01c9e7331acd49923f4bec3e91e57892a7ef2685/PanelMethod%200.2.0.py
[0.1.0]: https://github.com/AidenOgle/2DPanelMethod/blob/c8168b5cd490e402faf819ab6691faab28032685/PanelMethod%200.1.0.py
