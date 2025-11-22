import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation
from matplotlib import colors
from modules.bodysetup import*                  # module containing body setup functions
from modules.vortexprofile import*              # module containing Vortex Blob method functions
from modules.plotting import*                   # module containing plotting functions



# Define Simulation Parameters
chord = 1
density = 1000                          # kg/m3
nu = 1e-6                               # kinematic viscosity, m^2/s
num_panels = 50                         # body subdivisions
simlength = 5                           # length of sim in seconds, floats supported
sample_freq = 100                       # sampling frequency
condition = 'kelvin'                    # closure condition (only 'kelvin' is accurately supported, kutta conditions depreciated)
wakepoint_offset = [0, 0]               # offset distance for new wakepoints from TE to avoid singularities 
reference_frame = 'body'                # 'fluid' | 'body' , which frame is stationary
coeff_pressure_method = 'expanded'      # 'expanded' | 'simple' , allows for expanded form of Cp definition which avoids subtractive cancellation for small induced velocity values


# Body Parameters
body_params = {
    'bodytype': 'NACA',                 # 'line | 'ellipse' | 'circle' | 'NACA' | 'custom'
    'customx': None,                    # custom x values, 1 x num_panels+1 array
    'customy': None,                    # custom y values, 1 x num_panels+1 array
    'NACA': '0012',                     # e.g.:     '0012' | '4415' | '6412' | etc.
    'quarterchord': 0.25*chord,         # Set to zero to define rotation about the body center. Only functional for 'line' 'ellipse', & 'NACA'
    'eccentricity': 0.95,               # 0 to 1 inclusive, only functional for 'ellipse'
    'trailingedge': 'closed',           # Correction factor for trailing edge thickness. 'closed' | 'open' | 'prescribed'
    'bodyshape': 'closed'               # 'closed' or 'open' bodyshape. Used in Kutta condition enforcement. Required for 'custom' bodyshape, automatically assigned for others
}

# Pitching Motion Parameters
pitching_params = {
    'pitchingtype': 'constant',         # 'constant' | 'ramp' | 'ramp_rate' | 'periodic' | 'custom'
    'a1': 5,                            # initial angle of attack (deg)
    'a2': 15,                           # final angle of attack (deg)           (only used in 'between' method)
    'rate': None,                       # angluar pitch rate (deg/s)            (only used in 'rate' method)
    'amplitude': 15,                    # amplitude of periodic pitching (deg)  (only used in 'periodic' method)
    'period': simlength/2,                # pitching period (s)                   (only used in 'periodic' method)   
    'phase': 0                          # pitching phase offset (deg)           (only used in 'periodic' method) 
}


# Vortex Core Regularization Parameters
blobmethod = 'kaufmann'                 # 'none' | 'kaufmann' | 'lamb-oseen' | 'vatistas'  
vortexcore_kernel = 0.01                          # adjustable coefficient to fine-tune core size scale
spreading = False
# scale = np.mean(panel_lengths)          # core size length scale (reccomend either 'np.mean(panel_lengths)' or 'chord')        
# core_initial = scale * kernel

# Wake Vortex Decay Parameters
decaymethod = 'none'                                    # 'none' | 'exponential' | 'diffusive' (requires vortex spreading)
decay_kernel = 0.5                                      # adjustable coefficient for decay reference distance (larger -> slower decay)
tau = simlength * decay_kernel                          # decay timescale in seconds. used for 'exponential' method
# r_eff = np.mean(panel_lengths) * decay_kernel           # effective radius. used for 'diffusive' method


# Flow Velocity Parameters
u_initial = 1
u_final = 1


# Nondimensionalize time and compute the number of steps based on sample_freq and simlength
u_ref = np.max([u_initial, u_final])
timestep_target = chord / (u_ref * sample_freq)         # t* = (c/V_free(0))/Hz
num_step = int(simlength/timestep_target)
delT = simlength/num_step


# Initialize body geometry
[xinitial_body, yinitial_body, num_panels, bodyshape] = body_creation(num_panels=num_panels, chord=chord, **body_params)        # Call statement for body_creation function

body_position_initial = np.column_stack((xinitial_body, yinitial_body))                     # shape: (num_panels+1, 2)
vortex_position_initial = points_along_panels(body_position_initial, 0.25)
colloc_position_initial = points_along_panels(body_position_initial, 0.75)


# Initialize pitching values at timesteps
aoa_step = pitching(num_step=num_step, delT=delT,**pitching_params)            # Call statement for pitching function


# Displacement due to Pitching
rotation_matrix = rotation(num_step, aoa_step)                            # Defines rotation matrix at every timestep based on pitching profile
body_positions = np.einsum('ij,njk->nik', body_position_initial, rotation_matrix)               # shape: (num_step, num_panels+1, 2)
vortex_positions = np.einsum('ij,njk->nik', vortex_position_initial, rotation_matrix)           # shape: (num_step, num_panels, 2)
colloc_positions = np.einsum('ij,njk->nik', colloc_position_initial, rotation_matrix)           # shape: (num_step, num_panels, 2)


# Define body normals
panel_vectors = np.diff(body_positions, axis=1)                             # shape: (num_step, num_panels, 2)
panel_lengths = np.linalg.norm(panel_vectors, axis=2)                       # shape: (num_step, num_panels)

panel_dx = np.abs(panel_vectors[:, :, 0])                                   # shape: (num_step, num_panels)
panel_dy = np.abs(panel_vectors[:, :, 1])

panel_tangents = panel_vectors / panel_lengths[:, :, None]                  # shape: (num_step, num_panels, 2)
body_normals = perp(panel_tangents)                                        # shape: (num_step, num_panels, 2)
body_normals = - body_normals


# Dependent Vortex Parameters
vortexcore_scale = np.mean(panel_lengths)          # core size length scale (reccomend either 'np.mean(panel_lengths)' or 'chord')        
core_initial = vortexcore_scale * vortexcore_kernel
r_eff = np.mean(panel_lengths) * decay_kernel           # effective radius. used for 'diffusive' method

if blobmethod == 'none':
    spreading = False
    core_initial = 0

core2_body = core_initial**2


# Construct Freestream Velocity
freestream_velocity = np.zeros((num_step, 2))
freestream_velocity[:, 0] = np.linspace(u_initial, u_final, num_step)
freestream_velocity[:, 1] = 0
freestream_field = np.tile(freestream_velocity[:, None, :], (1, num_panels, 1))


# Displacement due to flow. Used if switching from a stationary body to a stationary fluid frame later
flow_displacement = np.cumsum(freestream_velocity, axis=0) * delT


# Preallocation for cached arrays
wake_positions = np.zeros((num_step, num_step, 2))
core2_wake = np.zeros((num_step, num_step)) 
gamma_body = np.zeros((num_step, num_panels + 1))           # vortex strength of bound body vortices and the newly shed wake vortex
gamma_wake_original = np.zeros((num_step, num_step))
gamma_wake_effective = np.zeros((num_step, num_step))
gamma_body_history = np.zeros((2, num_panels))             # single timestep cache for use in 'unsteady' calculation

induced_velocity_field_body = np.zeros((num_step, num_panels, 2)) 
velocity_field_body = np.zeros((num_step, num_panels, 2)) 
velocity_field_wake = np.zeros((num_step, num_step, 2))
coeff_pressure = np.zeros((num_step, num_panels))
coeff_lift = np.zeros((num_step))
coeff_lift_kj = np.zeros((num_step))
coeff_drag = np.zeros((num_step))  



# Process Timer initialization
start_time = time.process_time()

# Initialize optional progress bar
display_progress = True                 # toggles runtime progress bar display in terminal. Slight performance cost but helpful in debugging    
bar_length = 50
idler_interval = 0.1                    # time interval (seconds) btwn idle animation frame updates
idlerloop = ['⠁ ','⠃ ','⠇ ','⠧ ','⠷ ','⠿ ','⠿⠁','⠿⠃','⠿⠇','⠿⠧', '⠿⠷', '⠿⠿']
    # Wall-time idle animation: visual indicator of activity / smooth progress. tied to process clock instead of timestep
        # Alternative animations
        # ['◜','◝','◞','◟']
        # ['|', '/', '‒', '—', '‒', '\\']
        # ['⠁ ','⠃ ','⠇ ','⠧ ','⠷ ','⠿ ','⠿⠁','⠿⠃','⠿⠇','⠿⠧', '⠿⠷', '⠿⠿']
        # ['▁','▂','▃','▄','▅','▆','▇','█','▇','▆','▅','▄','▃','▂','▁'] 
display_process_time = True

#============================================================================================================================================================================================

# Main Timestep Loop
for i in range(num_step):

    num_wakepoints = i + 1

    # Access precomputed values for current timestep
    body = body_positions[i, :, :]              # (num_panels+1, 2)
    body_vortex = vortex_positions[i, :, :]     # (num_panels, 2)
    colloc = colloc_positions[i, :, :]          # (num_panels, 2)
    normal = body_normals[i, :, :]              # (num_panels, 2)


    # Wakepoint location initialization and reindexing
    if i > 0:
        # Orders indexing of wake values as a queue (newest=0 to oldest=-1)
        gamma_wake_original[i, 1:i+1] = gamma_wake_original[i - 1, 0:i]
        gamma_wake_effective[i, 1:i+1] = gamma_wake_effective[i - 1, 0:i]
        wake_positions[i, 1:i+1, :] = wake_positions[i - 1, 0:i, :]


    # Initialize newly shed wakepoint position at the TE
    wake_positions[i, 0, :] = body[-1, :] + wakepoint_offset


    # Preallocation for A & B Matrices
    A = np.zeros((num_panels+1, num_panels+1))
    B = np.zeros((num_panels+1))

    # Update wake vortex core size for 'spreading' case
    wakepoint_age = delT * np.arange(num_wakepoints)
    if spreading == True:
        core2_wake[i, :num_wakepoints] = core_initial**2 + 4 * nu * wakepoint_age
    else:
        core2_wake[i, :num_wakepoints] = core_initial**2


    # induced velocity coeff matrix at colloc points due to body vortices
    body_body_induced_vel_coeff = inducedvelocity(colloc, 
                                                  body_vortex, 
                                                  core2_body, 
                                                  blobmethod)
    
    # induced velocity coeff matrix at colloc points due to wake vortices
    body_wake_induced_vel_coeff = inducedvelocity(colloc, 
                                                  wake_positions[i, :num_wakepoints, :], 
                                                  core2_wake[i, :num_wakepoints], 
                                                  blobmethod)

    # Assemble system of equations
        # scaled_normals = normal * panel_lengths[i, :, None]         # shape: (num_panels, 2)
    A[:-1, :-1] = np.einsum('ijk,ik->ij', body_body_induced_vel_coeff, normal)          # body ind vel contribution
    A[:-1, -1] = np.einsum('ij,ij->i', body_wake_induced_vel_coeff[:, 0, :], normal)    # new wakepoint ind vel contribution

    B[:-1] = -np.einsum('ij,j->i', normal, freestream_velocity[i])                      # freestream vel contribution
    if i > 0:                                                                           # previous wake ind vel contribution
        wake_ind_vel = body_wake_induced_vel_coeff[:, 1:, :] * gamma_wake_effective[i, 1:num_wakepoints][None, :, None]
        B[:-1] += -np.sum(np.einsum('ijk,ik->ij', wake_ind_vel, normal), axis=1)                                                                              

    # Apply closure condition for system of equations
    match condition:
        case 'kelvin':
            # Kelvin's Circulation Theorem
            A[-1, :] = 1   
            if i == 0:                             
                B[-1] = 0  
            else:                                     
                B[-1] = -np.sum(gamma_wake_original[i, 1:num_wakepoints])   
        case 'kutta_thin':
            # Kutta Condition for single surface config
            A[-1, -2] = 1
            A[-1, -1] = 1
            B[-1] = 0
        case 'kutta_thick_unsteady':
            # Kutta Condition for dual surface config, wake-TE balancing
            A[-1, 0] = 1
            A[-1, -2] = 1 
            A[-1, -1] = 1
            B[-1] = 0 
        case 'kutta_thick_quasisteady':
            # Kutta Condition for dual surface config, TE balancing
            A[-1, 0] = 1
            A[-1, -2] = 1
            B[-1] = 0
    
    # Solve for unkown vortex strengths and store values appropriately
    gamma_body[i, :] = np.linalg.solve(A, B)
    gamma_body_history[1, :] = gamma_body_history[0, :]           # previous timestep gamma_body
    gamma_body_history[0, :] = gamma_body[i, :-1]                 # current timestep gamma_body
    gamma_wake_original[i, 0] = gamma_body[i, -1]                 # store new wakepoint strength
    gamma_wake_effective[i, 0] = gamma_body[i, -1]

    # Wake Vortex Decay 
    match decaymethod:
        case 'none':
            decayfactor = 1
        case 'exponential':
            decayfactor = np.exp(-delT/tau)
        case 'diffusive':
            if spreading == True:
                decayfactor = 1 - np.exp(-r_eff**2/core2_wake[i, 1:num_wakepoints])
            else:
                raise TypeError('Diffusive wake decay requires Wake Vortex Spreading to be enabled\n')
        case _:
            raise TypeError("Unrecognized Wake Vortex Decay Method. Supported methods are 'none', 'exponential', & 'diffusive'\n")
    
    gamma_wake_effective[i, 1:num_wakepoints] *= decayfactor                # update previously shed vortices strength

    # Find resultant velocity field at the collocation points
    induced_velocity_field_body[i] = (np.einsum('ijk,j->ik', body_body_induced_vel_coeff, gamma_body[i, :-1]) +                         # body induced vel contribution
                                      np.einsum('ijk,j->ik', body_wake_induced_vel_coeff, gamma_wake_effective[i, :num_wakepoints])     # wake induced vel contribution
                                      )
                            

    velocity_field_body[i] = freestream_field[i] + induced_velocity_field_body[i]       # total velocity field at body


    # Find local Pressure Coefficients Cp
    u_inf2 = np.sum(freestream_field[i]**2, axis=1)                     # redundant, repeated values for when freestream is spatially uniform. Left general for flexiblity
    # u_inf2 = np.sum(freestream_velocity[i]**2)
    match coeff_pressure_method:
        # Mathematically equivalent
        case 'simple':
            u_total2 = np.sum(velocity_field_body[i]**2,  axis=1)  
            coeff_pressure[i] = 1 - u_total2/(u_inf2 + 1e-15)               # 1e-15 to prevent division by zero
        case 'expanded':
            # Slightly more stable for very small values of induced velocity (u_ind ~1e-8)
            u_induced2 = np.sum(induced_velocity_field_body[i]**2, axis=1)
            free_dot_ind = np.sum(freestream_field[i] * induced_velocity_field_body[i], axis=1)
            coeff_pressure[i] = ((-2 * free_dot_ind) - u_induced2)/(u_inf2 + 1e-15)         # 1e-15 to prevent division by zero


    # Calculate steady Lift and Drag coeffs
    coeff_lift[i] = -np.sum(coeff_pressure[i, :] * normal[:, 1] * panel_lengths[i, :]) / chord
    coeff_drag[i] = np.sum(coeff_pressure[i, :] * normal[:, 0] * panel_lengths[i, :]) / chord  

    # Kutta-Joukowski Lift
    coeff_lift_kj[i] = -2*np.sum(gamma_body[i, :-1]) / (np.linalg.norm(freestream_velocity[i]) * chord)


    #==== Wake Roll-up ====#

    # induced velocity coeff matrix at wake points due to body vortices
    wake_body_induced_vel_coeff = inducedvelocity(wake_positions[i, :num_wakepoints, :],
                                              body_vortex,
                                              core2_body,
                                              blobmethod)

    # induced velocity coeff matrix at wake points due to wake vortices
    wake_wake_induced_vel_coeff = inducedvelocity(wake_positions[i, :num_wakepoints, :],
                                              wake_positions[i, :num_wakepoints, :],
                                              core2_wake[i, :num_wakepoints],
                                              blobmethod)

    # Eliminate wake-wake self-influence terms
    np.fill_diagonal(wake_wake_induced_vel_coeff[:, :, 0], 0.0)             # redundant for some core smoothing methods
    np.fill_diagonal(wake_wake_induced_vel_coeff[:, :, 1], 0.0)

    # Find resultant velocity field at the wake points
    velocity_field_wake[i, :num_wakepoints, :] = (freestream_velocity[i][None, :] +                                                                 # freestream contribution (broadcasting for efficiency)
                                                  np.einsum('ijk,j->ik', wake_body_induced_vel_coeff, gamma_body[i, :-1]) +                         # body ind vel contribution
                                                  np.einsum('ijk,j->ik', wake_wake_induced_vel_coeff, gamma_wake_effective[i, :num_wakepoints])     # wake ind vel contribution
                                                  )

    # Advect wakepoints using the total wake velocity field (induced + freestream)
    wake_positions[i, :num_wakepoints, :] += velocity_field_wake[i, :num_wakepoints, :]  * delT


    # Progress Bar
    if display_progress and (i + 1) %  max(1, num_step // 100) == 0:
        progress_fraction = (i + 1) / num_step
        filled_segments = int(round(progress_fraction * bar_length))
        progress_percent = f"{int(progress_fraction * 100)} %"
        progress_bar = "█" * filled_segments + "░" * (bar_length - filled_segments)

        idle_elapsed = time.time() - start_time
        idler_index = int(idle_elapsed / idler_interval) % len(idlerloop)
        idler_frame = idlerloop[idler_index]

        print(f"\r{idler_frame} Progress: |{progress_bar}| {progress_percent} ", end="", flush=True)

#============================================================================================================================================================================================
# Diagnostics


#===============================================================================

# Switch reference frames
if reference_frame == 'fluid':
    arrays = [body_positions, vortex_positions, colloc_positions, wake_positions, body_normals]
    for arr in arrays:
        arr -= flow_displacement[:, np.newaxis, :]      # Make sure shapes broadcast correctly


# Toggle visibility of different elements in the animated plot
display_body = True
display_bound = True
display_colloc = False
display_wake = True
display_normal = False
display_spreading = False

auto_xaxis_scaling = True

# Toggle body forces plots
plot_forces = True


# Dynamic scaling of X axis
if auto_xaxis_scaling:
    xaxis_min = np.min(body_positions[:, :, 0]) - 0.25*chord
    xaxis_max = np.max(wake_positions[:, :, 0]) + 0.25*chord
elif reference_frame == 'fluid':
    xaxis_min = -5.5
    xaxis_max = 1
elif reference_frame == 'body':
    xaxis_min = -0.75
    xaxis_max = 5.75
main_width = xaxis_max - xaxis_min

# Y axis scaling
yaxis_min = -1              # Current configuration has y-axis prescribed. Could be made dynamic using same method as x-axis
yaxis_max = 1
main_height = yaxis_max - yaxis_min

# Figure Sizing Parameters
plotscalefactor = 1.5       # Used to scale plot to be a reasonable size
colorbar_width = 0.2
padding = 0.5
fig_width = main_width*plotscalefactor + colorbar_width + padding
fig_height = main_height*plotscalefactor + colorbar_width + padding

# Plot Setup and Gridspace layout
fig = plt.figure(1, figsize=(fig_width, fig_height))
gs = gridspec.GridSpec(2, 2, figure=fig, width_ratios=[main_width*plotscalefactor, colorbar_width], 
                       height_ratios=[main_height*plotscalefactor, colorbar_width], hspace=0.2, wspace=0.05)
ax = fig.add_subplot(gs[0, 0])
ax.set_aspect('equal', adjustable='box', anchor='SW')
ax.set_title("Bound & Wake Vortex Evolution", fontsize=14, pad=10)

cax_right = fig.add_subplot(gs[0, 1])
cax_under = fig.add_subplot(gs[1, 0])

# Normalize colorbar values against vortex strengths
norm_bound = colors.CenteredNorm(vcenter=0, halfrange=np.percentile(np.abs(gamma_body), 99))
norm_wake = colors.CenteredNorm(vcenter=0, halfrange=np.percentile(np.abs(gamma_wake_effective), 99))
    # colorbars normalized independently, clips 1% outliers
cmap_bound = 'RdBu_r'       # colormap for vortex strengths
cmap_wake = 'RdBu_r'        # can choose separate colormaps for clarity (ex: 'PiYG')


# default marker sizes
bound_markersize = 5                      
wakepoint_markersize = np.full((num_step, num_step), 10.0)

# Initialize plot elements
bodyplot, = ax.plot([], [], 'k', lw=1 ,animated=True, zorder=1)
xdata_body, ydata_body = [], []

collocplot, = ax.plot([], [], 'c.', animated=True, zorder=2)
xdata_colloc, ydata_colloc = [], []

boundplot = ax.scatter([], [], s=bound_markersize, c=[], cmap=cmap_bound, animated=True, norm=norm_bound, zorder=3)
xdata_bound, ydata_bound = [], []

wakeplot = ax.scatter([], [], s=wakepoint_markersize[0, 0], c=[], cmap=cmap_wake, animated=True, norm=norm_wake, zorder=4)
xdata_wake, ydata_wake = [], []

normalplot = ax.quiver(colloc_position_initial[:,0], colloc_position_initial[:,1], np.full(num_panels, np.nan), np.full(num_panels, np.nan), color="r", scale=25, width=0.001)
xdata_normal, ydata_normal = [], []

time_box = ax.text(0.02, 0.95, '', transform=ax.transAxes, fontsize=8, fontfamily='Consolas', verticalalignment='center', horizontalalignment='left')
velocity_box = ax.text(0.02, 0.90, '', transform=ax.transAxes, fontsize=8, fontfamily='Consolas', verticalalignment='center', horizontalalignment='left')
aoa_box = ax.text(0.02, 0.85, '', transform=ax.transAxes, fontsize=8, fontfamily='Consolas', verticalalignment='center', horizontalalignment='left')


# Update markersize for 'spreading'
spreading_displayscale = 0.01                                           # scaling factor for percent change in marker size
if spreading and display_spreading:
    wakepoint_markersize *= ((core2_wake - core_initial**2) / (core_initial**2) * spreading_displayscale + 1)        # percent change with scaling factor
    wakepoint_markersize = np.clip(wakepoint_markersize, 0, None)       # sanitize array
    # Nonphysical size, shows scaled percent difference in coresize against a prescribed initial markersize
    # cant use actual coresize for markersize because order of ~1e-9 to ~1e-5 -> very small with large relative difference -> hard to meaningfully scale


# Create and format colorbars in the figure
cbar_norm_mode = 'mixed'        # 'mixed' | 'boundscaled' | 'wakescaled'
update_colorbars(
    fig=fig,
    right_ax=cax_right,
    under_ax=cax_under,
    boundplot=boundplot,
    wakeplot=wakeplot,
    norm_bound=norm_bound,
    norm_wake=norm_wake,
    mode=cbar_norm_mode          
)


# Plot Initialization
def init():
    ax.set_xlim(xaxis_min, xaxis_max)
    ax.set_ylim(yaxis_min, yaxis_max)
    time_box.set_text('')
    velocity_box.set_text('')
    aoa_box.set_text('')
    return [bodyplot, collocplot, boundplot, wakeplot, normalplot, time_box, velocity_box, aoa_box]

# Main plot updating function
def update(frame):

    plotdisplay = []

    # Body Plotting
    if display_body:
        bodyplot.set_data(body_positions[frame, :, 0], body_positions[frame, :, 1])
        plotdisplay.append(bodyplot)
    
    # Collocation Point Plotting
    if display_colloc:
        collocplot.set_data(colloc_positions[frame, :, 0], colloc_positions[frame, :, 0])
        plotdisplay.append(collocplot)

    # Bound Vortex Point Plotting
    if display_bound:
        boundplot.set_offsets(vortex_positions[frame])
        boundplot.set_array(gamma_body[frame])
        plotdisplay.append(boundplot)

    # Wake Vortex Point Plotting
    if display_wake:
        wkpts = frame + 1
        wakeplot.set_offsets(wake_positions[frame, :wkpts, :])
        wakeplot.set_array(gamma_wake_effective[frame, :wkpts])
        wakeplot.set_sizes(wakepoint_markersize[frame, :wkpts])
        # wakeplot.set_norm(colors.CenteredNorm())
        # colorbar.update_normal(wakeplot)
        plotdisplay.append(wakeplot)

    # Panel Normal Plotting
    if display_normal:
        normalplot.set_offsets(colloc_positions[frame])
        normalplot.set_UVC(body_normals[frame, :, 0], body_normals[frame, :, 1])
        plotdisplay.append(normalplot)

    # Display current simulation time (nondimensional)
    time_box.set_text(f"t = {frame * delT:.2f}")
    plotdisplay.append(time_box)

    # Display current (freestream) velocity
    u_current = freestream_velocity[frame, 0]
    velocity_box.set_text(f"U = {u_current:.2f}")
    plotdisplay.append(velocity_box)

    # Display current AOA
    aoa_current = -aoa_step[frame]
    aoa_box.set_text(f"α = {aoa_current:.2f}")
    plotdisplay.append(aoa_box)


    return plotdisplay



# Animation Display
ani = FuncAnimation(fig, update, frames=num_step, init_func=init, blit=True, interval=10, repeat_delay=1000)
    # 10ms between frames, 1000ms pause between playthroughs


# Scale body colorbar with main plot height
ax_pos = ax.get_position()
cb_bound_pos = cax_right.get_position()
cax_right.set_position([cb_bound_pos.x0, ax_pos.y0, cb_bound_pos.width, ax_pos.height])

# Scale wake colorbar with main plot width
cb_wake_pos = cax_under.get_position()
cax_under.set_position([ax_pos.x0, cb_wake_pos.y0, ax_pos.width, cb_wake_pos.height])

# secondary figure scaling after aspect ratio fixing and colorbar adjustments
gridtransform(gs, 0.95, 0, 0, [ax, cax_right, cax_under])        


# Body Forces Plot
if plot_forces:
    plt.figure(2)
    plt.title('Body Forces vs. Time')
    plotting_timestep = np.arange(num_step)[:] * delT
    plt.scatter(plotting_timestep, coeff_lift[:], s=5, label='Lift Coeff')
    plt.scatter(plotting_timestep, coeff_lift_kj[:], s=5, label='Lift Coeff (KJ)')
    plt.scatter(plotting_timestep, coeff_drag[:], s=5, label='Drag Coeff')
    aoa = -aoa_step[:]
    if body_params['bodytype'] == 'NACA' and body_params['NACA'][0] == '0':
        # Plots expected lift_coeff from thin airfoil theory
        # Limited to symmetric airfoils, cambered has an offset that varies by profile
        cl_expected = 2*np.pi * aoa
        plt.plot(plotting_timestep, cl_expected, c='k', label='Thin Airfoil C\u2097')
    plt.legend()
    plt.figure(3)
    plt.title('Body Forces vs. AOA')
    
    plt.scatter(aoa, coeff_lift[:], s=5, label='Lift Coeff')
    plt.scatter(aoa, coeff_lift_kj[:], s=5, label='Lift Coeff (KJ)')
    plt.scatter(aoa, coeff_drag[:], s=5, label='Drag Coeff')
    if body_params['bodytype'] == 'NACA' and body_params['NACA'][0] == '0':
        # Plots expected lift_coeff from thin airfoil theory
        # Limited to symmetric airfoils, cambered has an offset that varies by profile
        cl_expected = 2*np.pi * aoa
        plt.plot(aoa, cl_expected, c='k', label='Thin Airfoil C\u2097')
    plt.legend()

    plt.figure(1)

# Process Timer output
if display_process_time:
    end_time = time.process_time()
    elapsed = end_time - start_time
    print("\r" + " " * 100 + "\r", end="", flush=True)
    print(f"Process time: {elapsed:.4f} s\n")


# Toggle saving the animated plot
saveanimation = False
saveformat = 'gif'              # format
savepath = "./plots"            # location file is saved to
if saveanimation:
    animation_saving(ani, saveformat, "wakeplot", path=savepath)

# Toggle plotting all figures
display_plot = True
if display_plot:
    plt.show()



