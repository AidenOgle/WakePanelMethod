import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation
from matplotlib import colors
from modules.bodysetup import*              # module containing body setup functions
from modules.vortexblob import*             # module containing Vortex Blob method functions

# Process Timer initialization
start_time = time.process_time()

# Define Simulation Parameters
chord = 1
density = 1000                          # kg/m3
nu = 1e-6                               # kinematic viscosity, m^2/s
num_panels = 30
simlength = 5                           # length of sim in seconds, floats supported
subdivision = 100
num_step = int(simlength * subdivision)           # hundreths subdivision


# Flow Velocity Parameters
u_initial = -1
u_final = -1
u_step = np.linspace(u_initial, u_final, num_step)
delT = chord / (abs(u_step[0]) * subdivision)


# Body Parameters
body_params = {
    'bodytype': 'NACA',                 # 'line | 'ellipse' | 'circle' | 'NACA' | 'custom'
    'num_panels': num_panels,               # variable passing, do not modify
    'chord': chord,                         # variable passing, do not modify
    'customx': None,                    # custom x values, 1 x num_panels+1 array
    'customy': None,                    # custom y values, 1 x num_panels+1 array
    'NACA': '6412',                     # e.g.:     '0012' | '4415' | '6412'
    'quarterchord': 0.25*chord,         # Set to zero to define rotation about the body center. Only functional for 'line' 'ellipse', & 'NACA'
    'eccentricity': 0.95,               # 0 to 1 inclusive, only functional for 'ellipse'
    'trailingedge': 'closed',           # Correction factor for trailing edge thickness. 'closed' | 'open' | 'prescribed'
    'bodyshape': 'open'                 # 'closed' or 'open' bodyshape. Used in Kutta condition enforcement. Required for 'custom' bodyshape, automatically assigned for others
}

[xinitial_body, yinitial_body, num_panels, bodyshape] = body_creation(**body_params)        # Call statement for body_creation function

# Panel Length Components
panel_dx = np.abs(np.diff(xinitial_body))
panel_dy = np.abs(np.diff(yinitial_body))
panel_lengths = np.sqrt(panel_dx**2 + panel_dy**2)


# Pitching Motion Parameters
pitching_params = {
    'pitchingtype': 'periodic',         # 'constant' | 'between' | 'rate' | 'periodic' | 'custom'
    'num_step': num_step,                   # variable passing, do not modify
    'a1': 0,                            # initial angle of attack (deg)
    'a2': None,                         # final angle of attack (deg)           (only used in 'between' method)
    'rate': None,                       # angluar pitch rate (deg/s)            (only used in 'rate' method)
    'amplitude': 15,                    # amplitude of periodic pitching (deg)  (only used in 'periodic' method)
    'period': simlength,                # pitching period (s)                   (only used in 'periodic' method)   
    'phase': 0                          # pitching phase (deg)                  (only used in 'periodic' method) 
}

alpha_step = pitching(**pitching_params)            # Call statement for pitching function


# Vortex Blob Regularization Parameters
vortexblob_param = {
    'method': 'linear',                 # 'linear' | 'gaussian' | 'none'
    'scale': np.mean(panel_lengths),    # core size length scale (reccomend either 'np.mean(panel_lengths)' or 'const. * chord')        
    'spreading': False,                 # toggle for wake vortex spreading behavior
    'kernel': 0.5                       # adjustable coefficient to fine-tune core size
}

[blobmethod, spreading, core_initial] = vortexblobinit(**vortexblob_param)      # Call statement for Vortex Blob Method initialization


# Wake Vortex Decay Parameters
decaymethod = 'none'                                    # 'none' | 'exponential' | 'diffusive' (requires vortex spreading)
decay_kernel = 1                                        # adjustable coefficient for decay reference distance (larger -> slower decay)
tau = simlength * decay_kernel                          # decay timescale in seconds. used for 'exponential' method
r_eff = np.mean(panel_lengths) * decay_kernel           # effective radius. used for 'diffusive' method


# Displacement due to flow Precalculation
u_disp = np.zeros((num_step, 1))
for i in range(1, num_step):
    u_disp[i] = u_disp[i-1] + u_step[i-1]*delT      
    # nondimensionalized against initial flow


# Displacement due to pitching Precalculation
body_xrot = np.zeros((num_step, num_panels+1))
body_yrot = np.zeros((num_step, num_panels+1))
for i in range(num_step):
    body_xrot[i] = xinitial_body * np.cos(alpha_step[i]) - yinitial_body * np.sin(alpha_step[i])   
    body_yrot[i] = xinitial_body * np.sin(alpha_step[i]) + yinitial_body * np.cos(alpha_step[i])


# Preallocation for cached arrays
x_body, y_body = [np.zeros((num_step, num_panels + 1)) for _ in range(2)]                           # body geometry arrays (num_step, num_panels+1)
x_vortex, y_vortex, x_colloc, y_colloc = [np.zeros((num_step, num_panels)) for _ in range(4)]       # body fluid arrays (num_step, num_panels)
x_wake, y_wake, core2, gamma_wake = [np.zeros((num_step, num_step)) for _ in range(4)]              # wake point arrays (num_step, num_step) = (num_wakepoints, num_wakepoints)[final]
gamma_panel = np.zeros((num_step, num_panels + 1))
lift, drag = [np.zeros(num_step) for _ in range(2)]
coeff_lift, coeff_drag = [np.zeros(num_step) for _ in range(2)]

gamma_panel_history = np.zeros((2, num_panels))             # single timestep cache for use in 'unsteady' calculation

progress = "["+50*"_"+"]     "                              # Setup for progress bar

#============================================================================================================================================================================================

# Main Timestep Loop
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

    # Update location of trailing edge vorticies. Index wake info from newest [0] to oldest [i]
    if i > 0:
        gamma_wake[i, 1:i+1] = gamma_wake[i-1, 0:i]
        x_wake[i, 1:i+1] = x_wake[i-1, 0:i]
        y_wake[i, 1:i+1] = y_wake[i-1, 0:i]
    x_wake[i, 0] = (x_body[i, -1] - x_colloc[i, -1])                                                         \
        /np.sqrt((x_colloc[i, -1] - x_body[i, -1])**2 + (y_colloc[i, -1] - y_body[i, -1])**2)               \
        * delT * abs(u_step[i]) + x_body[i, -1]
    y_wake[i, 0] = (y_body[i, -1] - y_colloc[i, -1])                                                         \
        /np.sqrt((x_colloc[i, -1] - x_body[i, -1])**2 + (y_colloc[i, -1] - y_body[i, -1])**2)               \
        * delT * abs(u_step[i]) + y_body[i, -1]

    # Preallocation for Induced Velocity Coefficients and A & B Matrices
    num_wakepoints = i+1
    induced_vel_coeff = np.ones((num_panels, num_panels, 2))
    wake_induced_vel_coeff = np.ones((num_panels, num_wakepoints, 2))
    A = np.zeros((num_panels+1, num_panels+1))
    B = np.zeros((num_panels+1))

    # Calculate Pairwise Influence of Body-Body and Wake-Body points
    xc = x_colloc[i, :] 
    yc = y_colloc[i, :]
    xv = x_vortex[i, :]        
    yv = y_vortex[i, :]
    dx = xc.reshape(-1, 1) - xv.reshape(1, -1)
    dy = yc.reshape(-1, 1) - yv.reshape(1, -1)
    r2 = dx**2 + dy**2

    xw = x_wake[i, :num_wakepoints]
    yw = y_wake[i, :num_wakepoints]
    dxw = xc.reshape(-1, 1) - xw.reshape(1, -1)
    dyw = yc.reshape(-1, 1) - yw.reshape(1, -1) 
    r2_wake = dxw**2 + dyw**2

    # Wake Vortex Spreading
    wakepoint_age = delT * np.arange(num_wakepoints)
    if spreading == True:
        core2[i, :num_wakepoints] = core_initial**2 + 4 * nu * wakepoint_age
    else:
        core2[i, :num_wakepoints] = core_initial**2

    # Compute Body-Body Induced Velocity Coefficient Matrix
    induced_vel_coeff[:, :, 0] = dy / (2 * np.pi * r2)              # induced horizontal velocity coeff
    induced_vel_coeff[:, :, 1] = -dx / (2 * np.pi * r2)             # induced vertical velocity coeff

    # Vector Blob Method Regularization for Wake-Body Interaction
    blobfactor = vortexblobfactor(r2_wake, core2, i, blobmethod)

    # Compute Wake-Body Induced Velocity Coefficient Matrix
    wake_induced_vel_coeff[:, :, 0] = dyw * blobfactor              # induced horizontal velocity coeff
    wake_induced_vel_coeff[:, :, 1] = -dxw * blobfactor             # induced vertical velocity coeff
     
    # Build A and B matrices
    nx = normal_body[0, :]
    ny = normal_body[1, :]
    A[:-1, :-1] = induced_vel_coeff[:, :, 0] * nx.reshape(-1, 1) + induced_vel_coeff[:, :, 1] * ny.reshape(-1, 1)
    A[:-1, -1] = wake_induced_vel_coeff[:, 0, 0] * normal_body[0, :] + wake_induced_vel_coeff[:, 0, 1] * normal_body[1, :]
   
    B[:-1] = motion_body[0, :] * normal_body[0, :] + motion_body[1, :] * normal_body[1, :]      # Nonpenetration condition
    B[:-1] -= np.sum(wake_induced_vel_coeff[:, 1:, 0] * gamma_wake[i, 1:num_wakepoints] * nx.reshape(-1, 1) + wake_induced_vel_coeff[:, 1:, 1] * gamma_wake[i, 1:num_wakepoints] * ny.reshape(-1, 1), axis=1) 
   
    # Enforce Kutta Condition
    match bodyshape:
        case 'open':
            A[-1, :] = 1                                         
            B[-1] = -np.sum(gamma_wake[i, 1:num_wakepoints])        # Global Circulation constraint (weak form)
        case 'closed':
            A[-1, 0] = -1
            A[-1, -1] = 1
            B[-1] = 0                                               # Strong Form Kutta condition
    
    # Solve Circulation and update gamma terms
    gamma_panel[i, :] = np.linalg.solve(A, B)
    gamma_panel_history[1, :] = gamma_panel_history[0, :]           # previous timestep gamma_panel
    gamma_panel_history[0, :] = gamma_panel[i, :-1]                 # current timestep gamma_panel
    gamma_wake[i, 0] = gamma_panel[i, -1]                           # shed wake gamma

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
    gamma_wake[i, 1:num_wakepoints] *= decayfactor             # update previously shed vortices strength

    # Solve Lift, Drag, Coeff Lift, and Coeff Drag
    steady_lift = density * np.sum(gamma_panel_history[0, :]) * np.cos(-alpha_step[i]) * abs(u_step[i])
    steady_drag = density * np.sum(gamma_panel_history[0, :]) * -np.sin(-alpha_step[i]) * abs(u_step[i])

    diff_gamma_p = (gamma_panel_history[0, :] - gamma_panel_history[1, :]) / delT            
    unsteady_normal = density * np.sum(diff_gamma_p * panel_dx)
    unsteady_axial = density * np.sum(diff_gamma_p * panel_dy)
    unsteady_lift = unsteady_normal * np.cos(-alpha_step[i]) - unsteady_axial * np.sin(-alpha_step[i])
    unsteady_drag = unsteady_normal * np.sin(-alpha_step[i]) + unsteady_axial * np.cos(-alpha_step[i])

    lift[i] = steady_lift + unsteady_lift
    coeff_lift[i] = lift[i] / (0.5*density*chord*u_step[i]**2)
    drag[i] = steady_drag + unsteady_drag                        
    coeff_drag[i] = drag[i] / (0.5*density*chord*u_step[i]**2)
    
    # Wake Roll-up
    rollup_matrix = np.zeros((num_wakepoints, num_wakepoints + num_panels, 2))
    body_rollup_matrix = np.zeros((num_wakepoints, num_panels, 2))
    wake_rollup_matrix = np.zeros((num_wakepoints, num_wakepoints, 2))
    gamma_total = np.concatenate((np.transpose(gamma_panel[i, :num_panels]), gamma_wake[i, :num_wakepoints]))

    dx_rollup = xw.reshape(-1, 1) - xv.reshape(1, -1)
    dy_rollup = yw.reshape(-1, 1) - yv.reshape(1, -1)
    r2_rollup = dx_rollup**2 + dy_rollup**2

    body_rollup_matrix[:, :, 0] = dy_rollup / (2*np.pi * r2_rollup)
    body_rollup_matrix[:, :, 1] = -dx_rollup / (2*np.pi * r2_rollup)

    dxw_rollup = xw.reshape(-1, 1) - xw.reshape(1, -1)
    dyw_rollup = yw.reshape(-1, 1) - yw.reshape(1, -1)
    r2w_rollup = dxw_rollup**2 + dyw_rollup**2
    np.fill_diagonal(r2w_rollup, np.inf)            # eliminate self-influence terms

    # Vector Blob Method Regularization for Wake-Wake Interaction
    blobfactor_rollup = vortexblobfactor(r2w_rollup, core2, i, blobmethod)

    wake_rollup_matrix[:, :, 0] = dyw_rollup * blobfactor_rollup
    wake_rollup_matrix[:, :, 1] = -dxw_rollup * blobfactor_rollup

    rollup_matrix = np.concatenate((body_rollup_matrix, wake_rollup_matrix), axis=1)
    
    x_wake[i, :num_wakepoints] += np.sum(rollup_matrix[:, :, 0] * gamma_total, axis=1) * delT
    y_wake[i, :num_wakepoints] += np.sum(rollup_matrix[:, :, 1] * gamma_total, axis=1) * delT

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
colorbar_width = 0.2
padding = 0.5

fig_width = main_width*plotscalefactor + colorbar_width + 2*padding
fig_height = main_height*plotscalefactor + 2*padding

# Plot Setup
fig = plt.figure(figsize=(fig_width, fig_height))
gs = gridspec.GridSpec(1, 2, width_ratios=[main_width*plotscalefactor, colorbar_width])
ax = fig.add_subplot(gs[0, 0])
ax.set_aspect('equal', adjustable='box')
cax = fig.add_subplot(gs[0, 1])
color_max = max(np.max(np.abs(gamma_wake)), np.max(np.abs(gamma_panel)))

vortex_markersize = 5                      # default marker sizes
wakepoint_markersize = 10

bodyplot, = ax.plot([], [], 'k', animated=True, zorder=1)
xdata_body, ydata_body = [], []

collocplot, = ax.plot([], [], 'c.', animated=True, zorder=2)
xdata_colloc, ydata_colloc = [], []

vortexplot = ax.scatter([], [], s=vortex_markersize, c=[], cmap='RdBu_r', animated=True, vmax=color_max, vmin=-color_max, zorder=3)
xdata_vortex, ydata_vortex = [], []

wakeplot = ax.scatter([], [], s=wakepoint_markersize, c=[], cmap='RdBu_r', animated=True, vmax=color_max, vmin=-color_max, zorder=4)
xdata_wake, ydata_wake = [], []

spreading_displayscale = 500                                        # Set to Zero to disable Spreading Display
wakepoint_markersize = ((core2 - core_initial**2) / (core_initial**2) * spreading_displayscale + 1) * wakepoint_markersize    # percent change with scaling factor
wakepoint_markersize = np.clip(wakepoint_markersize, 0, None)       # sanitize array

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
    vortexplot.set_array(gamma_panel[frame])

    # Wake Point Plotting
    xdata_wake = np.flip(x_wake, 1)[frame, -frame:]
    ydata_wake = np.flip(y_wake, 1)[frame, -frame:]
    data = np.stack([xdata_wake, ydata_wake]).T
    wakeplot.set_offsets(data)
    wakeplot.set_array(np.flip(gamma_wake, 1)[frame, -frame:])
    wakeplot.set_sizes(np.flip(wakepoint_markersize, 1)[frame, -frame:])
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

gs.tight_layout(fig)

# Scale size of colorbar
cb_pos = cax.get_position()
ax_pos = ax.get_position()
cax.set_position([cb_pos.x0, ax_pos.y0, cb_pos.width, ax_pos.height])


plt.show()

