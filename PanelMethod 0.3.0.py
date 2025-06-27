import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation
from matplotlib import colors

# Process Timer initialization
start_time = time.process_time()

# Define Panel Parameters
chord = 1
density = 1000                          # kg/m3
nu = 1e-6                               # kinematic viscosity, m^2/s
num_panels = 30
simlength = 2                           # length of sim in seconds, floats supported
num_step = int(simlength*100)           # hundreths subdivision, truncates
#num_bodies = 1


# 4-Digit NACA airfoil generator. Used in 'body_creation' function for 'NACA' method
def generateNACAfoil(NACA, quarterchord=0.25*chord, trailingedge='closed'):
    """4-Digit NACA airfoil generator
    """
    # http://airfoiltools.com/airfoil/naca4digit  
    m = float(str(NACA)[0])/100                     # Max camber 
    p= float(str(NACA)[1])/10                       # Location of max camber
    t = float(str(NACA)[2:])/100                    # Max thickness as fraction of chord
    tc = [0.2969, -0.126, -0.3516, 0.2843, -0.1036]     # constants for thickness profile
    global num_panels

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
    return [xNACA, yNACA]


# Body Definition
def body_creation(bodytype, customx=None, customy=None, NACA='0012', quarterchord=0.25*chord, eccentricity=0.95, trailingedge='closed', bodyshape='open'):
    """Defines geometry of body
    
    bodytype: 'line, 'ellipse', 'circle', 'NACA' or 'custom'
    customx: 1 x num_panels+1 array
    customy: 1 x num_panels+1 array
    NACA: 4-Digit NACA airfoil designation. Inputed as either an int or string
    quarterchord: set to zero to define rotation about the body center. Only functional for 'line' 'ellipse', & 'NACA'
    eccentricity: 0 to 1 inclusive. only functional for 'ellipse'. defaults to 0.95
    trailingedge: correction factor for if a closed trailing edge is required for airfoil. Accepts arguments 'closed', 'open', & 'prescribed'
    """
    global num_panels, chord
    try:
        match bodytype:
            case 'line':
                xinitial_body = (np.linspace(-0.5*chord, 0.5*chord, num_panels + 1) + quarterchord)
                yinitial_body = np.zeros((num_panels + 1))
                bodyshape = 'open'

            case 'ellipse':
                param = np.linspace(0, 2*np.pi, num_panels + 1)
                if not 0 <= eccentricity <= 1:
                    raise ValueError("Invalid eccentricity. expected value between 0 and 1")
                a = 0.5*chord
                b = a*np.sqrt(1-eccentricity**2)
                xinitial_body = a*np.cos(param) + quarterchord
                yinitial_body = b*np.sin(param)
                bodyshape = 'closed'

            case 'circle':
                param = np.linspace(0, 2*np.pi, num_panels + 1)
                xinitial_body = np.cos(param) * 0.5*chord
                yinitial_body = np.sin(param) * 0.5*chord
                bodyshape = 'closed'

            case "NACA":
                [xinitial_body, yinitial_body] = generateNACAfoil(NACA, quarterchord, trailingedge)
                bodyshape = 'closed'

            case 'custom':
                if np.shape(customx) == (num_panels+1, ) and np.shape(customy) == (num_panels+1, ):
                    xinitial_body = customx
                    yinitial_body = customy

                else:
                    print("\nImproper custom body array size(s)\n")
            case _:
                print("\nImproper body type selected\n")
        return [xinitial_body, yinitial_body, bodyshape]
    except TypeError:
        print("\nArguments missing or improperly entered in 'body_creation' function:\n")


[xinitial_body, yinitial_body, bodyshape] = body_creation('NACA', NACA=6412)    # Call statement for body_creation function

"""
Examples of different arguments for 'body_creation':
    body_creation('line', quarterchord=0)
    body_creation('ellipse', eccentricity=0.80)
    body_creation('NACA', NACA=6412)          # cambered airfoil
    body_creation('NACA', NACA='0015')        # symmetric airfoil. NACA number must be formatted as a string
    body_creation('NACA', NACA=4415, trailingedge='prescribed', quarterchord=0)
    body_creation('custom', np.linspace(0, 1, num_panels+1), np.zeros(num_panels+1))
"""

# Flow Velocity Parameters Definition
u_initial = -1
u_final = -1
u_step = np.linspace(u_initial, u_final, num_step)
delT = 0.01*chord/abs(u_step[0])


# Displacement due to flow
u_disp = np.zeros((num_step, 1))
for i in range(1, num_step):
    u_disp[i] = u_disp[i-1] + u_step[i-1]*delT      
    # nondimensionalized against initial flow


# Pitching Parameters Definition
def pitching(pitchingtype, a1=0, a2=None, rate=None, amplitude=None, period=simlength, phase=0):
    """Defines pitch angle for each time step.
    
    constant    ->  pitching('constant', a1=angle*)  
    between     ->  pitching('between', a1=angle*, a2=angle)             
    rate        ->  pitching('rate', a1=angle*, rate=ratevalue)                
    periodic    ->  pitching('periodic', a1=angle*, amplitude=value, period=value, phase=value*)
    custom      ->  pitching('custom', [1 x num_step array])
        - *optional parameters, 'a1' & 'phase' defaults to 0  to have horizontal start with rising leading edge
        - 'angle', 'rate', 'amplitude', & 'phase' values are in degrees, 'period' is in nondimensionalized seconds
    """
    try:
        match pitchingtype:
            case 'constant':
                alpha_initial = -np.radians(a1)
                return alpha_initial * np.ones((num_step, 1))
            
            case 'between':
                alpha_initial = -np.radians(a1)
                alpha_final = -np.radians(a2)
                return np.linspace(alpha_initial, alpha_final, num_step)
            
            case 'rate':
                alpha_initial = -np.radians(a1)
                omega = -np.radians(rate)*num_step/100       # find final position based on timescale
                return np.linspace(alpha_initial, omega, num_step)
            
            case 'periodic':
                alpha_initial = -np.radians(a1)
                amp = np.radians(amplitude)
                per = 2*np.pi/period/100
                phaseshift = -np.radians(phase)                           
                return alpha_initial + amp*np.sin(-per*np.arange(num_step) + phaseshift)
            
            case 'custom':
                if np.shape(a1) == (num_step, ):
                    return -np.radians(a1)
                else:
                    print("\nImproper custom pitching profile array size\n")
            case _:
                print("\nImproper pitching type selected\n")
    except TypeError:
        print("\nArguments missing or improperly entered in 'pitching' function:\n")

    
alpha_step = pitching('periodic', a1=0, amplitude=15, period=simlength)           # Call statement for pitching function

"""
Examples of different arguments for 'pitching':
    pitching('periodic', a1=15, amplitude=45, period=5, phase=180)
    pitching('rate', rate=30)
    pitching('constant')
    pitching('custom', np.linspace(0, 90, num_step))
    pitching('periodic', a1=0, amplitude=15, period=simlength)
"""

# Displacement due to pitching
body_xrot = np.zeros((num_step, num_panels+1))
body_yrot = np.zeros((num_step, num_panels+1))
for i in range(num_step):
    body_xrot[i] = xinitial_body * np.cos(alpha_step[i]) - yinitial_body * np.sin(alpha_step[i])   
    body_yrot[i] = xinitial_body * np.sin(alpha_step[i]) + yinitial_body * np.cos(alpha_step[i])


# Find Panel Length Components
panel_dx = np.abs(np.diff(xinitial_body))
panel_dy = np.abs(np.diff(yinitial_body))
panel_lengths = np.sqrt(panel_dx**2 + panel_dy**2)
#x_mid = 0.5 * (xinitial_body[:-1] + xinitial_body[1:])
#chordwise_order = np.argsort(x_mid)         # gives panel indices in order form LE to TE

# Vortex Blob Regularization Method
blobmethod = 'none'                 # 'linear' | 'gaussian' | 'none'
spreading = False                   # wake vortex spreading
blob_kernel = 0.1                   # tuneable coefficient for vortex core size
core_initial = np.mean(panel_lengths) * blob_kernel

# Wake Vortex Decay
decaymethod = 'none'                # 'exponential' | 'diffusive' (requires vortex spreading) | 'none'
decay_kernel = 0.5                  # tuneable coefficient for decay reference distance (larger -> slower decay)
tau = simlength * decay_kernel      # decay timescale in seconds
r_eff = np.mean(panel_lengths) * decay_kernel


# Preallocation for cached arrays
x_body, y_body = [np.zeros((num_step, num_panels + 1)) for _ in range(2)]                       # body geometry arrays (num_step, num_panels+1)
x_vortex, y_vortex, x_colloc, y_colloc = [np.zeros((num_step, num_panels)) for _ in range(4)]   # body fluid arrays (num_step, num_panels)
x_wake, y_wake, core2, gamma_w = [np.zeros((num_step, num_step)) for _ in range(4)]             # wake point arrays (num_step, num_step) = (num_wakepoints, num_wakepoints)[final]
circulation = np.zeros((num_step, num_panels + 1))
lift, drag = [np.zeros(num_step) for _ in range(2)]
coeff_lift, coeff_drag = [np.zeros(num_step) for _ in range(2)]

gamma_p = np.zeros((2, num_panels))              # single timestep cache for use in 'unsteady' calculation

progress = "["+50*"_"+"]     "                  # Setup for progress bar

for i in range(num_step):
    # Total displacement of the body
    x_body[i, :] = body_xrot[i] + u_disp[i]
    y_body[i, :] = body_yrot[i]

    # Build vortex points
    x_vortex[i, :] = (x_body[i, 1:] - x_body[i, 0:-1]) * 0.25 + x_body[i, 0:-1]
    y_vortex[i, :] = (y_body[i, 1:] - y_body[i, 0:-1]) * 0.25 + y_body[i, 0:-1]

    # Build collocation points
    x_colloc[i, :] = (x_body[i, 1:] - x_body[i, 0:-1]) * 0.75 + x_body[i, 0:-1]
    y_colloc[i, :] = (y_body[i, 1:] - y_body[i, 0:-1]) * 0.75 + y_body[i, 0:-1]

    # Build vectors normal to the body
    normal_body = np.zeros((2, num_panels))
    normal_body[0, :] = (y_colloc[i, :] - y_vortex[i, :]) / np.sqrt((y_colloc[i, :] - y_vortex[i, :])**2 + (x_colloc[i, :] - x_vortex[i, :])**2)
    normal_body[1, :] = -(x_colloc[i, :] - x_vortex[i, :]) / np.sqrt((y_colloc[i, :] - y_vortex[i, :])**2 + (x_colloc[i, :] - x_vortex[i, :])**2)

    # Build vectors for body motion
    motion_body = np.zeros((2, num_panels))
    motion_body[0, :] = u_step[i]
    motion_body[1, :] = 0

    # Update location of trailing edge vorticies
    if i > 0:
        gamma_w[i, 1:i+1] = gamma_w[i-1, 0:i]
        x_wake[i, 1:i+1] = x_wake[i-1, 0:i]
        y_wake[i, 1:i+1] = y_wake[i-1, 0:i]
    x_wake[i, 0] = (x_body[i, -1] - x_colloc[i, -1])                                                         \
        /np.sqrt((x_colloc[i, -1] - x_body[i, -1])**2 + (y_colloc[i, -1] - y_body[i, -1])**2)               \
        * delT * abs(u_step[i]) + x_body[i, -1]
    y_wake[i, 0] = (y_body[i, -1] - y_colloc[i, -1])                                                         \
        /np.sqrt((x_colloc[i, -1] - x_body[i, -1])**2 + (y_colloc[i, -1] - y_body[i, -1])**2)               \
        * delT * abs(u_step[i]) + y_body[i, -1]

    # Preallocation for Induced Velocity Coefficients and A & B Matrices (Biot-Savart Law calculations)
    num_wakepoints = i+1
    induced_vel_coeff = np.ones((num_panels, num_panels, 2))
    wake_induced_vel_coeff = np.ones((num_panels, num_wakepoints, 2))
    A = np.zeros((num_panels+1, num_panels+1))
    B = np.zeros((num_panels+1))

    # Calculate Pairwise Influence of Body-Body and Wake-Body points
    xc = x_colloc[i, :].reshape(-1, 1)        # shape (num_panels, 1)
    yc = y_colloc[i, :].reshape(-1, 1)
    xv = x_vortex[i, :].reshape(1, -1)        # shape (1, num_panels)
    yv = y_vortex[i, :].reshape(1, -1)
    dx = xc - xv
    dy = yc - yv
    r2 = dx**2 + dy**2

    xw = x_wake[i, :num_wakepoints].reshape(1, -1)
    yw = y_wake[i, :num_wakepoints].reshape(1, -1)
    dxw = xc - xw
    dyw = yc - yw 
    r2_wake = dxw**2 + dyw**2

    wakepoint_age = delT * np.arange(num_wakepoints)

    # Wake Vortex Spreading
    if spreading == True:
        core2[i, :num_wakepoints] = core_initial**2 + 4 * nu * wakepoint_age
    else:
        core2[i, :num_wakepoints] = core_initial**2

    # Vector Blob Method Regularization
    match blobmethod:
        case 'none':
            blobfactor = 1 / r2_wake                   
        case 'linear':
            blobfactor = 1 / (r2_wake + core2[i, :num_wakepoints])
        case 'gaussian':
            blobfactor = (1 - np.exp(-r2_wake/core2[i, :num_wakepoints])) / r2_wake
        case _:
            raise TypeError("Unrecognized Vortex Blob Method. Supported methods are 'none', 'linear', & 'gaussian'\n")

    # Compute Body Induced Velocity Coefficient Matrix
    induced_vel_coeff[:, :, 0] = dy / (2 * np.pi * r2)
    induced_vel_coeff[:, :, 1] = -dx / (2 * np.pi * r2)

    # Compute Wake Induced Velocity Coefficient Matrix
    wake_induced_vel_coeff[:, :, 0] = dyw / (2*np.pi) * blobfactor
    wake_induced_vel_coeff[:, :, 1] = -dxw / (2*np.pi) * blobfactor
     
    # Build A and B matrices
    nx = normal_body[0, :].reshape(-1, 1)
    ny = normal_body[1, :].reshape(-1, 1)
    A[:-1, :-1] = induced_vel_coeff[:, :, 0] * nx + induced_vel_coeff[:, :, 1] * ny
    A[:-1, -1] = wake_induced_vel_coeff[:, 0, 0] * normal_body[0, :] + wake_induced_vel_coeff[:, 0, 1] * normal_body[1, :]
   
    B[:-1] = motion_body[0, :] * normal_body[0, :] + motion_body[1, :] * normal_body[1, :]      # Nonpenetration condition
    B[:-1] -= np.sum(wake_induced_vel_coeff[:, 1:, 0] * gamma_w[i, 1:num_wakepoints] * nx + wake_induced_vel_coeff[:, 1:, 1] * gamma_w[i, 1:num_wakepoints] * ny, axis=1) 
   
    # Enforce Kutta Condition
    match bodyshape:
        case 'open':
            A[-1, :] = 1                                         
            B[-1] = -np.sum(gamma_w[i, 1:num_wakepoints])       # Global Circulation constraint (weak form)
        case 'closed':
            A[-1, 0] = -1
            A[-1, -1] = 1
            B[-1] = 0                                           # Strong Form Kutta condition
    
    # Solve Circulation and update gamma terms
    circulation[i, :] = np.linalg.solve(A, B)
    gamma_p[1, :] = gamma_p[0, :]                     # previous timestep gamma_p
    gamma_p[0, :] = circulation[i, :-1]               # current timestep gamma_p
    gamma_w[i, 0] = circulation[i, -1]                # shed wake vorticity

    # Wake Vortex Decay
    match decaymethod:
        case 'none':
            decayfactor = 1
        case 'exponential':
            decayfactor = np.exp(-delT/tau)
        case 'diffusive':
            if spreading == True:
                decayfactor = 1 - np.exp(-r_eff**2/core2[i, 1:num_wakepoints])
            else:
                raise TypeError('Diffusive wake decay requires Wake Vortex Spreading to be enabled\n')
        case _:
            raise TypeError("Unrecognized Wake Vortex Decay Method. Supported methods are 'none', 'exponential', & 'diffusive'\n")
    gamma_w[i, 1:num_wakepoints] *= decayfactor          # update previously shed vortices strength

    # Solve Lift, Drag, Coeff Lift, and Coeff Drag
    steady_lift = density * np.sum(gamma_p[0, :]) * np.cos(-alpha_step[i]) * abs(u_step[i])
    steady_drag = density * np.sum(gamma_p[0, :]) * -np.sin(-alpha_step[i]) * abs(u_step[i])

    diff_gamma_p = (gamma_p[0, :] - gamma_p[1, :]) / delT            
    unsteady_normal = density * np.sum(diff_gamma_p * panel_dx)
    unsteady_axial = density * np.sum(diff_gamma_p * panel_dy)
    unsteady_lift = unsteady_normal * np.cos(-alpha_step[i]) - unsteady_axial * np.sin(-alpha_step[i])
    unsteady_drag = unsteady_normal * np.sin(-alpha_step[i]) + unsteady_axial * np.cos(-alpha_step[i])

    lift[i] = steady_lift + unsteady_lift
    coeff_lift[i] = lift[i] / (0.5*density*chord*u_step[i]**2)
    drag[i] = steady_drag + unsteady_drag                        
    coeff_drag[i] = drag[i] / (0.5*density*chord*u_step[i]**2)

    # Wake Roll-up
    rollup_matrix = np.ones((num_wakepoints, num_wakepoints+num_panels, 2))
    x_total = np.append(x_vortex[i, :], x_wake[i, :num_wakepoints])
    y_total = np.append(y_vortex[i, :], y_wake[i, :num_wakepoints])
    gamma_total = np.append(np.transpose(circulation[i, :num_panels]), gamma_w[i, :num_wakepoints])

    for j in range(num_wakepoints):                        # for each wake point that currently exists
        for k in range(num_wakepoints+num_panels):             # compare to all other wake points and vortex points
            dist_rollup = np.sqrt((x_wake[i, j] - x_total[k])**2 + (y_wake[i, j] - y_total[k])**2 + (0.04*chord)**2)
            rollup_matrix[j, k, 0] = 1/2/np.pi / dist_rollup * (y_wake[i, j] - y_total[k]) / dist_rollup
            rollup_matrix[j, k, 1] = 1/2/np.pi / dist_rollup * (-x_wake[i, j] + x_total[k]) / dist_rollup
            
            """
            rollup_matrix[j, k, 0] = 1/2/np.pi / np.sqrt((x_wake[j] - x_total[k])**2 + (y_wake[j] - y_total[k])**2 + (0.04*chord)**2)      \
                  * (y_wake[j] - y_total[k]) / np.sqrt((x_wake[j] - x_total[k])**2 + (y_wake[j] - y_total[k])**2 + (0.04*chord)**2)
            rollup_matrix[j, k, 1] = 1/2/np.pi / np.sqrt((x_wake[j] - x_total[k])**2 + (y_wake[j] - y_total[k])**2 + (0.04*chord)**2)      \
                  * (-x_wake[j] + x_total[k]) / np.sqrt((x_wake[j] - x_total[k])**2 + (y_wake[j] - y_total[k])**2 + (0.04*chord)**2)
            """

        x_wake[i, j] = x_wake[i, j] + delT*sum(rollup_matrix[j, :, 0]*gamma_total)
        y_wake[i, j] = y_wake[i, j] + delT*sum(rollup_matrix[j, :, 1]*gamma_total)

    # Progress Bar
    if (i+1) % (num_step//50) == 0:
        print("Progress: ", progress, end="", flush=True)
        print("\r", end="", flush=True)
        percent = str(int(((i+1)/num_step)*100)) + " %"
        progress = progress[:(i+1)//(num_step//50)] + "█" + progress[(i+1)//(num_step//50)+1:-len(percent)] + percent


# Dynamic scaling of plot axes
auto_xaxis_scaling = True
if auto_xaxis_scaling == True:
    xaxis_min = np.min(x_body) - 0.25*chord
    xaxis_max = np.max(x_body) + 0.25*chord
else:
    xaxis_min = -5.5
    xaxis_max = 1
main_width = xaxis_max - xaxis_min

yaxis_min = -1              # Current configuration has y-axis prescribed. Could be made dynamic using same method as x-axis
yaxis_max = 1
main_height = yaxis_max - yaxis_min

# Figure Sizing Parameters
plotscalefactor = 1.5       # Used to scale plot to be a reasonable size
colorbar_width = 0.1
padding = 0.75
wspace_inch = 0.1           # spacing between plot and color bar
wspace = 2 * wspace_inch / (main_width + colorbar_width)

fig_width = 2*padding + main_width*plotscalefactor + wspace_inch + colorbar_width
fig_height = main_height*plotscalefactor

# Plot Setup
fig = plt.figure(figsize=(fig_width, fig_height))         
fig.subplots_adjust(left=padding/fig_width, right=1-padding/fig_width)
gs = gridspec.GridSpec(1, 2, width_ratios=[main_width, colorbar_width], wspace=wspace)
ax = fig.add_subplot(gs[0, 0])
ax.set_aspect('equal')
cax = fig.add_subplot(gs[0, 1])
color_max = max(np.max(np.abs(gamma_w)), np.max(np.abs(circulation)))    

bodyplot, = ax.plot([], [], 'k', animated=True, zorder=1)
xdata_body, ydata_body = [], []

collocplot, = ax.plot([], [], 'c.', animated=True, zorder=2)
xdata_colloc, ydata_colloc = [], []

vortexplot = ax.scatter([], [], s=5, c=[], cmap='RdBu_r', animated=True, vmax=color_max, vmin=-color_max, zorder=3)
xdata_vortex, ydata_vortex = [], []

wakeplot = ax.scatter([], [], s=10, c=[], cmap='RdBu_r', animated=True, vmax=color_max, vmin=-color_max, zorder=4)
xdata_wake, ydata_wake = [], []

colorbar = fig.colorbar(wakeplot, cax=cax)

# Toggle visibility of different elements in the animated plot
display_body = True
display_vortex = True
display_colloc = False
display_wake = True

# Plot Initialization
def init():
    ax.set_xlim(xaxis_min, xaxis_max)
    ax.set_ylim(-1, 1)
    return bodyplot, 

# Main plot updating function
def update(frame):
    # Body Plotting
    xdata_body = x_body[frame]
    ydata_body = y_body[frame]
    bodyplot.set_data(xdata_body, ydata_body)
    
    # Collocation Point Plotting
    xdata_colloc = x_colloc[frame]
    ydata_colloc = y_colloc[frame]
    collocplot.set_data(xdata_colloc, ydata_colloc)

    # Vortex Point Plotting
    xdata_vortex = x_vortex[frame]
    ydata_vortex = y_vortex[frame]
    data_vortex = np.stack([xdata_vortex, ydata_vortex]).T
    vortexplot.set_offsets(data_vortex)
    vortexplot.set_array(circulation[frame])

    # Wake Point Plotting
    xdata_wake = np.flip(x_wake, 1)[frame, -frame:]
    ydata_wake = np.flip(y_wake, 1)[frame, -frame:]
    data = np.stack([xdata_wake, ydata_wake]).T
    wakeplot.set_offsets(data)
    wakeplot.set_array(np.flip(gamma_w, 1)[frame, -frame:])
    wakeplot.set_norm(colors.CenteredNorm())
    colorbar.update_normal(wakeplot)

    plotdisplay = []
    if display_body == True:
        plotdisplay.append(bodyplot)
    if display_vortex == True:
        plotdisplay.append(vortexplot)
    if display_colloc == True:
        plotdisplay.append(collocplot)
    if display_wake == True:
        plotdisplay.append(wakeplot)
    return plotdisplay

# Process Timer output
end_time = time.process_time()
print("\r"+" "*100, end="", flush=True)
print("\r", end="", flush=True)
print("\nProcess time: ", end_time-start_time, " s\n")

# Animation Display
ani = FuncAnimation(fig, update, frames=num_step, init_func=init, blit=True, interval=10, repeat_delay=1000)
    # 10ms between frames, 1000ms pause between playthroughs

# Scale size of colorbar
cb_pos = cax.get_position()
ax_pos = ax.get_position()
cax.set_position([cb_pos.x0, ax_pos.y0, cb_pos.width, ax_pos.height])

plt.show()

