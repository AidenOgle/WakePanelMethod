"""vortexblob.py

This module defines functions for vortex blob method initialization and coefficient computation

Functions:
- vortexblobinit(...)
- vortexblobfactor(...)
"""

import numpy as np

# Vortex Blob Regularization Method
def vortexblobinit(method, scale, spreading, kernel):
    """ Used to initialize Vortex Blob Regularization Method and gatekeep vortex spreading mechanism.
    
    method:     'linear' | 'gaussian' | 'none'
    scale:      length scale vortex core is defined against
    spreading:  toggle for wake vortex spreading behavior
    kernel:     adjustable coefficient to control core size
    """
    match method:
        case 'none':
            spreading = False
            core_initial = 0
        case 'linear':
            core_initial = scale * kernel
        case 'gaussian':
            core_initial = scale * kernel
        case _:
            raise TypeError("Unrecognized Vortex Blob Method. Supported methods are 'none', 'linear', & 'gaussian'\n")
    return [method, spreading, core_initial]


# Vortex Blob Method coefficient
def vortexblobfactor(r2, coresize, i, blobmethod):
    """Defines the coefficients used in the calculation of induced velocity coefficient matrices."""
    num_wakepoints = i + 1
    c2 = coresize[i, :num_wakepoints]
    match blobmethod:
        case 'none':
            coeff = 1 / (2 * np.pi * r2)                   
        case 'linear':
            coeff = 1 / (2 * np.pi * (r2 + c2))
        case 'gaussian':
            coeff = (1 - np.exp(-r2 / c2)) / (2 * np.pi * r2)
        case _:
            raise TypeError("Unrecognized Vortex Blob Method. Supported methods are 'none', 'linear', & 'gaussian'\n")
    return coeff


