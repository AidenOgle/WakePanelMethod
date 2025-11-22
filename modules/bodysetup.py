"""bodysetup.py

This module defines initialization functions for body geometry definition and body motion

Functions:
- generateNACAfoil(...)
- body_creation(...)
- pitching(...)
"""

import numpy as np

# 4-Digit NACA airfoil generator. Used in 'body_creation' function for 'NACA' method
def generateNACAfoil(NACA, num_panels, chord, quarterchord, trailingedge):
    """4-Digit NACA airfoil generator
    """
    # http://airfoiltools.com/airfoil/naca4digit  
    m = float(str(NACA)[0])/100                     # Max camber 
    p= float(str(NACA)[1])/10                       # Location of max camber
    t = float(str(NACA)[2:])/100                    # Max thickness as fraction of chord
    tc = [0.2969, -0.126, -0.3516, 0.2843, -0.1036]     # constants for thickness profile


    if num_panels % 2 == 1:                         # Odd number of panels
        print("\nOnly even numbers of panels are supported for NACA method. Increasing 'num_panels' from ", num_panels, " to ", num_panels+1, "\n")
        num_panels = num_panels + 1
    
    num_NACA = (num_panels)//2 + 1
    param = np.linspace(0, np.pi, num_NACA)
    xc = 0.5*(1 + np.cos(param))                    # cosine spacing (refinement around leading/trailing edges)
  
    yc = np.zeros(num_NACA)   
    dycdx = np.zeros(num_NACA)  

    match trailingedge:
        case 'open':
            tc[4] = -0.1015
        case 'closed':                              # correction for closed trailing edge
            tc[4] = -0.1036                         
        case 'prescribed':                          
            pass                                    # handled later after yt is defined
        case _:
            raise ValueError("Invalid input for 'trailingedge'. Expected either 'open', 'closed', or 'prescribed'")
    
    yt = 5*t*(tc[0]*xc**0.5 + tc[1]*xc + tc[2]*xc**2 + tc[3]*xc**3 + tc[4]*xc**4)       # thickness profile
    
    if trailingedge == 'prescribed':
        yt[-1] = 0                                   # prescribes thickness to be 0 at trailing edge instead of calculating
    
    for i, xc[i] in enumerate(xc):
        if p == 0:                                  # symmetric airfoil case
            pass
        elif xc[i] <= p:
            yc[i] = m/p**2 * (2*p*xc[i] - xc[i]**2)                     # camber profile pre max camber
            dycdx[i] = (2*m/p**2)*(p-xc[i])
        else:
            yc[i] = m/(1-p)**2 * ((1-2*p) + 2*p*xc[i] - xc[i]**2)       # camber profile post max camber
            dycdx[i] = 2*m/(1-p)**2 * (p - xc[i])
    theta = np.arctan(dycdx)
    xu = xc - yt*np.sin(theta)          # upper profile
    yu = yc + yt*np.cos(theta)
    xl = xc + yt*np.sin(theta)          # lower profile
    yl = yc - yt*np.cos(theta)

    xNACA = (np.append(xu[0:-1], np.flip(xl)) - 0.5) * chord + quarterchord      # flip lower array to maintain counterclockwise definition
    yNACA = np.append(yu[0:-1], np.flip(yl)) * chord                      
    return [xNACA, yNACA, num_panels]


# Body Geometry Definition
def body_creation(bodytype, num_panels, chord, customx, customy, NACA, quarterchord, eccentricity, trailingedge, bodyshape):
    """Defines body geometry as discretized points in 2D space.
    
    bodytype:       'line | 'ellipse' | 'circle' | 'NACA' | 'custom'
    customx:        1 x num_panels+1 array
    customy:        1 x num_panels+1 array
    NACA:           4-Digit NACA airfoil designation. Inputed as either an int or string
    quarterchord:   Set to zero to define rotation about the body center. Only functional for 'line' 'ellipse', & 'NACA'
    eccentricity:   0 to 1 inclusive. only functional for 'ellipse'. defaults to 0.95
    trailingedge:   Correction factor for trailing edge thickness. Accepts arguments 'closed' | 'open' | 'prescribed'
    bodyshape:      'thick' or 'thin' bodyshape. Used in Kutta condition enforcement. Required for 'custom' bodyshape, automatically assigned for others
    """
    try:
        match bodytype:
            case 'line':
                xinitial_body = (np.linspace(-0.5*chord, 0.5*chord, num_panels + 1) + quarterchord)
                yinitial_body = np.zeros((num_panels + 1))
                bodyshape = 'thin'

            case 'ellipse':
                param = np.linspace(0, 2*np.pi, num_panels + 1)
                if not 0 <= eccentricity <= 1:
                    raise ValueError("Invalid eccentricity. expected value between 0 and 1")
                a = 0.5*chord
                b = a*np.sqrt(1-eccentricity**2)
                xinitial_body = a*np.cos(param) + quarterchord
                yinitial_body = b*np.sin(param)
                bodyshape = 'thick'

            case 'circle':
                param = np.linspace(0, 2*np.pi, num_panels + 1)
                xinitial_body = np.cos(param) * 0.5*chord
                yinitial_body = np.sin(param) * 0.5*chord
                bodyshape = 'thick'

            case "NACA":
                [xinitial_body, yinitial_body, num_panels] = generateNACAfoil(NACA, num_panels, chord, quarterchord, trailingedge)
                bodyshape = 'thick'

            case 'custom':
                if np.size(customx) == num_panels+1 and np.size(customy) == num_panels+1:
                    xinitial_body = customx
                    yinitial_body = customy
                else:
                    print("\nMismatch between number of points provided and number of panels. Adjusting num_panels to match number of points\n")
                    num_panels = np.size(customx) - 1
            case _:
                print("\nImproper body type selected\n")
        return [xinitial_body, yinitial_body, num_panels, bodyshape]
    except TypeError:
        print("\nArguments missing or improperly entered in 'body_creation' function:\n")


# Pitching Parameters Definition
def pitching(pitchingtype, num_step, delT, a1, a2, rate, amplitude, period, phase):
    """Defines pitch angle for each time step.
    
    constant    ->  pitching('constant', a1=angle)  
    ramp        ->  pitching('ramp', a1=angle, a2=angle)             
    ramp_rate   ->  pitching('ramp_rate', a1=angle, rate=ratevalue)                
    periodic    ->  pitching('periodic', a1=angle, amplitude=value, period=value, phase=value)
    custom      ->  pitching('custom', [1 x num_step array])
    """
    try:
        match pitchingtype:
            case 'constant':
                alpha_initial = -np.radians(a1)
                return alpha_initial * np.ones((num_step))
            
            case 'ramp':
                alpha_initial = -np.radians(a1)
                alpha_final = -np.radians(a2)
                return np.linspace(alpha_initial, alpha_final, num_step)
            
            case 'ramp_rate':
                alpha_initial = -np.radians(a1)
                omega = -np.radians(rate)*num_step/100       # find final position based on timescale
                return np.linspace(alpha_initial, omega, num_step)
            
            case 'periodic':
                alpha_initial = -np.radians(a1)
                amp = np.radians(amplitude)
                w = 2*np.pi/period
                phaseshift = -np.radians(phase)
                t = np.arange(num_step) * delT                           
                return alpha_initial + amp*np.sin(-w*t + phaseshift)
            
            case 'custom':
                if np.shape(a1) == (num_step, ):
                    return -np.radians(a1)
                else:
                    print("\nImproper custom pitching profile array size\n")
            case _:
                print("\nImproper pitching type selected\n")
    except TypeError:
        print("\nArguments missing or improperly entered in 'pitching' function:\n")


def rotation(num_step, alpha_step):
    """Rotation Transformation Matrix at every timestep
    """
    rotation = np.empty((num_step, 2, 2))                                       
    rotation[:, 0, 0] = np.cos(-alpha_step)
    rotation[:, 0, 1] = -np.sin(-alpha_step)
    rotation[:, 1, 0] = np.sin(-alpha_step)
    rotation[:, 1, 1] = np.cos(-alpha_step)
    return rotation


def points_along_panels(body, fraction):
    """Compute points along each panel at a given fraction

    Points are located at the fractional length along the panel in the downstream direction. 
    The freestream direction is important because collocation points are placed downstream of vortex points for numerical stability purposes
    Downstream is considered the positive x direction during definition.     
    """
    p0 = body[:-1, :]
    p1 = body[1:, :]
    contour = p1 - p0
    
    points = np.empty_like(p0)
    downstream_mask = (p1[:, 0] > p0[:, 0])         # checks if panel order is locally defined aligned with freestream
    points[downstream_mask] = p0[downstream_mask] + fraction * contour[downstream_mask]         # freestream aligned
    points[~downstream_mask] = p1[~downstream_mask] - fraction * contour[~downstream_mask]      # freestream counteraligned
    # points = p0 + fraction * contour
    return points
