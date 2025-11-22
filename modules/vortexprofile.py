"""vortexblob.py

This module defines functions for vortex blob method initialization and coefficient computation

Functions:
- perp(...)
- vortexblobinit(...)
- vortexblobfactor(...)
- inducedvelocity(...)
"""

import numpy as np

def perp(vecs):
    """Return 2D vectors rotated 90° CCW."""
    x, y = vecs[..., 0], vecs[..., 1]
    return np.stack([-y, x], axis=-1)

# # Vortex Blob Regularization Method
# def vortexblobinit(method, scale, spreading, kernel):
#     """ Used to initialize Vortex Blob Regularization Method and gatekeep vortex spreading mechanism.
    
#     method :    'linear' | 'gaussian' | 'none'
#     scale :     length scale vortex core is defined against
#     spreading : toggle for wake vortex spreading behavior
#     kernel :    adjustable coefficient to control core size
#     """
#     match method:
#         case 'none':
#             spreading = False
#             core_initial = 0
#         case 'kaufmann' | 'scully':
#             core_initial = scale * kernel
#         case 'lamb' | 'oseen' | 'lamb-oseen' | 'lamb oseen':
#             core_initial = scale * kernel
#         case _:
#             raise TypeError("Unrecognized Vortex Blob Method. Supported methods are 'none', 'linear', & 'gaussian'\n")
#     return [method, spreading, core_initial]


# Vortex Blob Method coefficient
def vortexblobfactor(r2, c2, blobmethod):
    """Regularizes vortex profiles used in induced velocity calculations to prevent singularities.
    Gives a geometric term that regularizes the vortex core profile for an effective radius.
    """
    match blobmethod:
        case 'none':
            coeff = 1 / (2 * np.pi * (r2 + 1e-12))
            # technically kaufmann, but with a very small fixed core size to prevent singularity                   
        case 'kaufmann' | 'scully':
            coeff = 1 / (2 * np.pi * (r2 + c2))
        case 'lamb' | 'oseen' | 'lamb-oseen' | 'lamb oseen':
            coeff = (1 - np.exp(-r2 / c2)) / (2 * np.pi * (r2 + 1e-12))
        case 'vatistas':
            n = 2
            coeff = 1 / (2 * np.pi * (r2**n + c2**n)**(1/n))
        case _:
            raise TypeError("Unrecognized Vortex Blob Method. Supported methods are 'none', 'linear', & 'gaussian'\n")
    return coeff


def inducedvelocity(probe_position, source_position, c2, blobmethod='kaufmann'):
    """Build geometric influence matrix v_ij for vortex induced velocity
    
    Parameters
    ----------
    probe : (n, 2) array
        Coordinates of probe points.
    source : (m, 2) array
        Coordinates of source points.
    c2 : (n, m) or similar
        vortex core size squared.
    blobmethod : str
        One of 'none', 'linear', or 'gaussian'.
    
    Returns
    -------
    v_ij : (n, m, 2) array
        Geometric influence coefficients.
    """
    diff = probe_position[:, None, :] - source_position[None, :, :]
    r2 = np.sum(diff**2, axis=-1)                                       # (n, m)
    blobfactor = vortexblobfactor(r2, c2, blobmethod)                   # (n, m)
    v_ij = perp(diff) * blobfactor[:, :, None]
    return v_ij




