import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

"""
2D Panel Method
Aiden Ogle, 2025
"""

# Define Panel Parameters
chord = 1
num_bodies = 1
num_panels = 25
num_step = 500      # hundreths of a second

# # Preallocate Arrays
# x_wake = np.ones((1,num_step))
# y_wake = np.ones((1,num_step))
# gamma_w = np.ones((1,num_step))
# gamma_p = np.zeros((2,num_panels))
# lift = np.ones((1,num_step))
# coeff_lift = np.ones((1,num_step))

# Body Definition
bodytype = 'ellipse'
quarterchord = 0.25*chord         # set to zero to define body about center
match bodytype:
    case 'line':
        xinitial_body = (np.linspace(-0.5*chord, 0.5*chord, num_panels + 1) + quarterchord)
        yinitial_body = np.zeros((num_panels + 1))  
    case 'ellipse':
        param = np.linspace(0,2*np.pi,num_panels + 1)
        eccentricity = 0.95
        a = 0.5*chord
        b = a*np.sqrt(1-eccentricity**2)
        xinitial_body = a*np.cos(param) + quarterchord
        yinitial_body = b*np.sin(param)
        """
        dist = np.sqrt((xinitial_body[1:]-xinitial_body[0:-1])**2 + (yinitial_body[1:]-yinitial_body[0:-1])**2)
        print(dist)
        
        # ellipse currently results in uneven panel legths. This is mostly noticable for low panel subdivisions
        # unclear how significant this will be, solution is nontrivial and would involve elliptic integral approximation
        """
    case 'circle':
        param = np.linspace(0,2*np.pi,num_panels + 1)
        xinitial_body = np.cos(param) * 0.5*chord
        yinitial_body = np.sin(param) * 0.5*chord



# Flow Velocity Parameters Definition
u_initial = -1
u_final = -1
u_step = np.linspace(u_initial,u_final,num_step)
delT = 0.01*chord/abs(u_step[0])

# Displacement due to flow
u_disp = np.zeros((num_step,1))
for i in range(1,num_step):
    u_disp[i] = u_disp[i-1] + u_step[i-1]*delT      
    # nondimensionalized against initial flow


# Pitching Parameters Definition
def pitching(args):
    """Defines prescribed pitch angle.
    
    constant pitch                      ->  pitching(angle)  
                                        ex: pitching(45)
    uniform pitch rate between values   ->  pitching([angle1, angle2])             
                                        ex: pitching([0, 90]) 
    uniform pitch rate using rate       ->  pitching(['rate', angle1, rate])                
                                        ex: pitching(['rate', 0, 90])
    periodic                            ->  pitching(['periodic',angle1,amplitude,period])
                                        ex: pitching(['periodic', 0, 45, 5])
    """
    if isinstance(args,(int, float)):
        alpha_initial = -np.radians(args)
        return alpha_initial * np.ones((num_step, 1))
    elif len(args) == 2 and isinstance(args[0],(int, float)):
        alpha_initial = -np.radians(args[0])
        alpha_final = -np.radians(args[1])
        return np.linspace(alpha_initial,alpha_final,num_step)
    elif args[0] == 'rate':
        alpha_initial = -np.radians(args[1])
        omega = -np.radians(args[2])*num_step/100
        return np.linspace(alpha_initial,omega,num_step)
    elif args[0] == 'periodic':
        alpha_initial = -np.radians(args[1])
        amplitude = np.radians(args[2])
        period = 2*np.pi/args[3]/100
        phase = np.radians(180)                           # method to change phase not currently implemented
        return alpha_initial + amplitude*np.sin(period*np.arange(num_step) + phase)
    else:
        print('Improper Pitching Parameters Input')
    # Likely will change how arguments are passed in this function to use *args or **kwargs


alpha_step = pitching(['periodic', 0, 45, 5])



# Displacement due to pitching
body_xrot = np.zeros((num_step,num_panels+1))
body_yrot = np.zeros((num_step,num_panels+1))
for i in range(num_step):
    body_xrot[i] = xinitial_body * np.cos(alpha_step[i]) - yinitial_body * np.sin(alpha_step[i])   
    body_yrot[i] = xinitial_body * np.sin(alpha_step[i]) + yinitial_body * np.cos(alpha_step[i])


x_body = np.zeros((num_step,num_panels+1))
y_body = np.zeros((num_step,num_panels+1))
x_vortex = np.zeros((num_step,num_panels))
y_vortex = np.zeros((num_step,num_panels))
x_colloc = np.zeros((num_step,num_panels))
y_colloc = np.zeros((num_step,num_panels))
normal_body = np.ones((num_step,2,num_panels))
motion_body = np.ones((num_step,2,num_panels))
for i in range(num_step):
    # Total displacement of the body
    x_body[i] = body_xrot[i] + u_disp[i]
    y_body[i] = body_yrot[i]

    # Build vortex points
    x_vortex[i] = (x_body[i,1:] - x_body[i,0:-1]) * 0.25 + x_body[i,0:-1]
    y_vortex[i] = (y_body[i,1:] - y_body[i,0:-1]) * 0.25 + y_body[i,0:-1]

    # Build collocation points
    x_colloc[i] = (x_body[i,1:] - x_body[i,0:-1]) * 0.75 + x_body[i,0:-1]
    y_colloc[i] = (y_body[i,1:] - y_body[i,0:-1]) * 0.75 + y_body[i,0:-1]

    # Build vectors normal to the body
    normal_body[i,0,:] = (y_colloc[i] - y_vortex[i])/np.sqrt((y_colloc[i] - y_vortex[i])**2 + (x_colloc[i] - x_vortex[i])**2)
    normal_body[i,1,:] = -(x_colloc[i] - x_vortex[i])/np.sqrt((y_colloc[i] - y_vortex[i])**2 + (x_colloc[i] - x_vortex[i])**2)

    # Build vectors for body motion
    motion_body[i,0,:] = u_step[i]
    motion_body[i,1,:] = 0



# Dynamic scaling of plot axes
auto_xaxis_scaling = True
if auto_xaxis_scaling == True:
    xaxis_min = np.min(x_body) - 0.25*chord
    xaxis_max = np.max(x_body) + 0.25*chord
else:
    xaxis_min = -5.5
    xaxis_max = 1
xaxis_size = xaxis_max - xaxis_min

yaxis_min = -1          # Current configuration has y-axis prescribed. Could be made dynamic using same method as x-axis
yaxis_max = 1
yaxis_size = yaxis_max - yaxis_min

axis_ratio = xaxis_size/yaxis_size
plotscalefactor = 3     # Used to scale plot to be a reasonable size. Lower if plot gets too big to prevent graphical warping

# Plot Setup
fig, ax = plt.subplots(figsize=(axis_ratio*plotscalefactor,plotscalefactor))
bodyplot, = ax.plot([], [], 'k')
xdata_body, ydata_body = [], []

vortexplot, = ax.plot([], [], 'r.')
xdata_vortex, ydata_vortex = [], []

collocplot, = ax.plot([], [], 'c.')
xdata_colloc, ydata_colloc = [], []

# Toggle visibility of different elements in the animated plot
display_body = True
display_vortex = True
display_colloc = False

# Plot Initialization
def init():
    ax.set_xlim(xaxis_min, xaxis_max)
    ax.set_ylim(-1, 1)
    return bodyplot,

# Main plot updating function
def update(frame):
    xdata_body = x_body[frame]
    ydata_body = y_body[frame]
    bodyplot.set_data(xdata_body, ydata_body)

    xdata_vortex = x_vortex[frame]
    ydata_vortex = y_vortex[frame]
    vortexplot.set_data(xdata_vortex, ydata_vortex)

    xdata_colloc = x_colloc[frame]
    ydata_colloc = y_colloc[frame]
    collocplot.set_data(xdata_colloc,ydata_colloc)

    plotdisplay = []
    if display_body == True:
        plotdisplay.append(bodyplot)
    if display_vortex == True:
        plotdisplay.append(vortexplot)
    if display_colloc == True:
        plotdisplay.append(collocplot)
    return plotdisplay


# Animation Display
ani = FuncAnimation(fig, update, frames=num_step, init_func=init, blit=True, interval=5, repeat_delay=1000)
    # 5ms between frames, 1000ms pause between playthroughs
plt.show()
