# CHANGELOG

All notable changes to this project will be documented in this file

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - YYYY/MM/DD

### Added

- Added a 4-Digit NACA airfoil option for body geometry. Supports both symmetric and cambered airfoils. Incorporates new argument 'trailingedge=' which allows the user to specifiy if the trailing edge has a finite thickness or is closed. Viable assignments are 'open', 'closed', and 'prescribed'. 'prescribed' and 'closed' should achieve the same thing, but 'closed' depends on an equation whereas 'prescribed' directly assigns the thickness at the trailing edge to be zero. This assignment may be depreciated in the future depending on usefulness

  - Due to the way vortex and collocation points are defined relative to body points based on order, the top and bottom surfaces have an offset in the order of vortex and collocation points compared to eachother. This is most noticeable for symmetric airfoils with a low number of panels (e.g. NACA='0015', num_panels=10, display_vortex=True, display_colloc=True)
  - Odd numbers of panels results in a duplicate body point at the leading edge (0,0), resulting in a coincident vortex and collocation point. The solution to this, if necessary, would require either a difference in number of panels between the upper and lower surfaces or a truncation of either the leading or tailing edge geometry

- Option to define custom body geometry using x-y coordinates. Currently inputted as two arrays, may change to have different syntax in the future to reflect how other programs export x-y data
- Option to define custom pitching profile using an array of values corresponding to timesteps. Allows for pitching profiles to be defined by functions/equations beyond the preconstructed ones provided

### Changed

- Reformatted the geometry creation section into a function structure as 'body_creation'
- Changed how arguments are passed to the 'pitching' function to be more consistent. Now uses a keyword argument structure
- Added an exception handling case for improper argument syntax and inputs for the 'pitching' and 'body_creation' functions to account for the user-unfriendliness of the keyword argument structure. Might implement prompt system in the future for easier user-input
- Added 'simlength' variable to drive 'num_step' to better reflect the time in whole units instead of hundreths. No direct usage yet, added primarily as a reference variable
  
### Fixed
  
## [0.1.0] -  2025/06/4
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

[0.1.0]: https://github.com/AidenOgle/2DPanelMethod/blob/c8168b5cd490e402faf819ab6691faab28032685/PanelMethod%200.1.0.py
