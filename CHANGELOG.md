# CHANGELOG

All notable changes to this project will be documented in this file

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - YYYY/MM/DD

## [0.2.0] - 2025/06/15

### Added

- FLUID SIMULATION: Incorporated vortex panel method used in original MATLAB code to find lift derived from body circulation and wake behavior. The main difference between the implementations in MATLAB and here is that the MATLAB implementation has local variables that are plotted and written over for each timestep iteration, whereas the implementation here stores values outside of the loop for each time step. This change leads to slightly decreased performance but better plotting stability and access to data at different time steps
  - BODY SIMULATION:
    - Have to double check if Lift (L) and the Coefficient of Lift (CL) are calculated correctly in accordance with nondimensionalization, as the original MATLAB version didn't allow for variable velocity or density
    - Intend to create animated plots for L and CL and generate data for other properties like drag
  - WAKE SIMULATION:
      - the MATLAB implementation structured the assignment of values for wake points as a queue, which works well for iterative plotting but makes cached plotting slightly obtuse. I intend to eventually change this structure to follow normal increasing array assignment for consistency, but this would require restructuring across multiple sections. Currently handled in the plotting section with a stopgap method of array reversal and indexing
      - Intend to implement prescribed viscocity effect for wake circulation dissipation
  - COMPLEXITY and OPTIMIZATION:
    - The Wake Roll-Up section forms the bulk of processing complexity being O(n<sup>3</sup>), with the body lift calculation being O(n<sup>2</sup>⋅m<sup>2</sup>) and the total complexity being O(n<sup>3</sup>+n<sup>2</sup>⋅m<sup>2</sup>), where 'n' is the number of time steps (and wake points) and 'm' is the number of panels.
    - Ideally, optimization would be accomplished through vectorization to eliminate extraneous loops. I do not know if this is possible with the current setup but I will explore this in the next version
    - Intend to implement parameters that allow wake calculation to be turned off or reduced. Reduction could be in the form of  time resolution by only generating/calculating wake points every nth time step, or by having a range of influence such that wake points are only compared against points that are close in index (i.e. +/-50 index)
  
- NACA AIRFOIL : Added a 4-Digit NACA airfoil option for body geometry to directly generate body profiles based on the 4-digit series airfoils. Supports both symmetric and cambered airfoils. [The equations on this site](http://airfoiltools.com/airfoil/naca4digit) were used as reference. This is an experimental feature, it may be more reliable to directly input x-y points using the 'custom' method

  - Incorporates new argument 'trailingedge=' which allows the user to specifiy if the trailing edge has a finite thickness or is closed. Viable assignments are 'open', 'closed', and 'prescribed'. 'prescribed' and 'closed' should achieve the same thing, but 'closed' uses an equation whereas 'prescribed' directly assigns the thickness at the trailing edge to be zero. This assignment may be depreciated in the future depending on usefulness
  - Due to the way vortex and collocation points are defined relative to body points based on order, the top and bottom surfaces have an offset in the order of vortex and collocation points compared to eachother. This is most noticeable for symmetric airfoils with a low number of panels (e.g. NACA='0015', num_panels=10, display_vortex=True, display_colloc=True)
  - Odd numbers of panels would cause a duplicate point at either the leading or trailing edges, resulting in coincident vortex and collocation points. To avoid this, the number of panels is shunted up to an even value. If an odd number of panels is desired, it is reccomended to use the 'custom' method  to manually define points

- CUSTOM BODY : Option to define custom body geometry using x-y coordinates. Currently inputted as two arrays, may change to have different syntax in the future to reflect how other programs export x-y data

- CUSTOM PITCHING : Option to define custom pitching profile using an array of values corresponding to timesteps. Allows for pitching profiles to be defined by functions/equations beyond the preconstructed ones provided

- Added a progress bar and program timer. Unlike the original MATLAB code, the results are not shown as they are generated and are instead displayed asynchronously after the program has finished. This setup improves the viewability of the results animation but leads to a longer wait time between initialization and output, especially for large panel numbers. The progress bar and program timer were implemented for quality of life and future optimization purposes

### Changed

- Reformatted the geometry creation section into a function structure as 'body_creation'
- Changed how arguments are passed to the 'pitching' function to be more consistent. Now uses a keyword argument structure
- Added an exception handling case for improper argument syntax and inputs for the 'pitching' and 'body_creation' functions to account for the user-unfriendliness of the keyword argument structure. Might implement prompt system in the future for easier user-input
- Added 'simlength' variable to drive 'num_step' to better reflect the time in whole units instead of hundreths. Allows floats. No direct usage yet, added primarily as a reference variable
  
### Fixed

- Fixed 'periodic' case for 'pitching' function to reflect direction conventions. The positive anglular direction for the body reference frame cooresponds to the negative global anglular direction from cartesian frame, which was not originally implemented correctly. Changed 'phase' argument to default to 0 to reflect intended default motion under currected system.

## [0.1.0] -  2025/06/04
Created initial version of Panel Method code
Only handles geometry, motion (displacement & pitching), and plotting animation so far. Wake shedding and pressure profile will be addressed in next update

### Added

- GEOMETRY:
  - created presets for 'line', 'ellipse', and 'circle' bodies
  - currently bodies are initially centered at the origin (either about their center or quarter-chord). Plan on allowing initial offsets in the future
- DISPLACEMENT:
  - changed from MATLAB code to allow for irregular displacements. Uses stepwise motion
  - currently only allows for only constant velocity or linear accel/deceleration
  - planned to change to follow same structure as PITCHING section and allow for more predefined or custom velocity profiles
- PITCHING:
  - added predefined pitching profiles 'constant', 'constant rate between values', 'constant rate using rate', and 'periodic'
  - currently structured as a function that passes 'args', will be changed in the future to be more consistent between cases
- PLOTTING:
  - uses FuncAnimation module from matplotlib library
  - added option to automatically scale x-axis to fit the displacement of the body
  - added plotscaling calculation to keep consistent scale between axes and prevent graphical distortion of the body
    - if the auto axis scaling causes the plotsize to exceed the display size, the scaling becomes inconsistent. Will fix next update
  -added options to display/hide the body, vortex points, and collocation points respectively

### Changed
  
### Fixed

[0.2.0]: 
[0.1.0]: https://github.com/AidenOgle/2DPanelMethod/blob/c8168b5cd490e402faf819ab6691faab28032685/PanelMethod%200.1.0.py
