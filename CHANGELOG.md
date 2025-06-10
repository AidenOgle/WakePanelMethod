# CHANGELOG

All notable changes to this project will be documented in this file

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - YYYY/MM/DD


### Added

### Changed
  
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
