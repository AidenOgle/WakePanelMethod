import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import colors

"""
2D Panel Method
Aiden Ogle, 2025
"""

start_time = time.process_time()

# Define Panel Parameters
chord = 1
density = 1000                  # kg/m3, water (not sure if this is what was intended in original matlab code, will doube check later)
num_bodies = 1
num_panels = 30
simlength = 2.5                         # length of sim in seconds, floats supported
num_step = int(simlength*100)           # hundreths subdivision


# 4-Digit NACA airfoil generator. Used in 'body_creation' function for 'NACA' method
def generateNACAfoil(NACA,quarterchord=0.25*chord,trailingedge='closed'):
    """4-Digit NACA airfoil generator
    """
    # http://airfoiltools.com/airfoil/naca4digit  
    m = float(str(NACA)[0])/100                     # Max camber 
    p= float(str(NACA)[1])/10                       # Location of max camber
    t = float(str(NACA)[2:])/100                    # Max thickness as fraction of chord
    tc = [0.2969,-0.126,-0.3516,0.2843,-0.1036]     # constants for thickness profile
    global num_panels

    if num_panels % 2 == 1:                         # Odd number of panels
        print("\nOnly even numbers of panels are supported for NACA method. Increasing 'num_panels' from ",num_panels," to ",num_panels+1,"\n")
        num_panels = num_panels + 1
    
    num_NACA = (num_panels)//2 + 1
    param = np.linspace(0,np.pi,num_NACA)
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
    
    for i,xc[i] in enumerate(xc):
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

    xNACA = (np.append(xu[0:-1],np.flip(xl)) - 0.5) * chord + quarterchord      # flip lower array to maintain counterclockwise definition
    yNACA = np.append(yu[0:-1],np.flip(yl)) * chord                      
    return [xNACA,yNACA]


# Body Definition
def body_creation(bodytype,customx=None,customy=None,NACA='0015',quarterchord=0.25*chord,eccentricity=0.95,trailingedge='closed'):
    """Defines geometry of body
    
    bodytype: 'line,'ellipse','circle', 'NACA' or 'custom'
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

            case 'ellipse':
                param = np.linspace(0,2*np.pi,num_panels + 1)
                if not 0 <= eccentricity <= 1:
                    raise ValueError("Invalid eccentricity. expected value between 0 and 1")
                a = 0.5*chord
                b = a*np.sqrt(1-eccentricity**2)
                xinitial_body = a*np.cos(param) + quarterchord
                yinitial_body = b*np.sin(param)

            case 'circle':
                param = np.linspace(0,2*np.pi,num_panels + 1)
                xinitial_body = np.cos(param) * 0.5*chord
                yinitial_body = np.sin(param) * 0.5*chord

            case "NACA":
                [xinitial_body,yinitial_body] = generateNACAfoil(NACA,quarterchord,trailingedge)

            case 'custom':
                if np.shape(customx) == (num_panels+1,) and np.shape(customy) == (num_panels+1,):
                    xinitial_body = customx
                    yinitial_body = customy

                else:
                    print("\nImproper custom body array size(s)\n")
            case _:
                print("\nImproper body type selected\n")
        return [xinitial_body,yinitial_body]
    except TypeError:
        print("\nArguments missing or improperly entered in 'body_creation' function:\n")


[xinitial_body,yinitial_body] = body_creation('NACA',NACA='6412')    # Call statement for body_creation function

"""
Examples of different arguments for 'body_creation':
    body_creation('line',quarterchord=0)
    body_creation('ellipse',eccentricity=0.80)
    body_creation('NACA',NACA=6412)          # cambered airfoil
    body_creation('NACA',NACA='0015')        # symmetric airfoil
    body_creation('NACA',NACA=4415,trailingedge='prescribed',quarterchord=0)
    body_creation('custom',np.linspace(0,1,num_panels+1),np.zeros(num_panels+1))
"""

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
def pitching(pitchingtype,a1=0,a2=None,rate=None,amplitude=None,period=None,phase=0):
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
                return np.linspace(alpha_initial,alpha_final,num_step)
            
            case 'rate':
                alpha_initial = -np.radians(a1)
                omega = -np.radians(rate)*num_step/100       # find final position based on timescale
                return np.linspace(alpha_initial,omega,num_step)
            
            case 'periodic':
                alpha_initial = -np.radians(a1)
                amp = np.radians(amplitude)
                per = 2*np.pi/period/100
                phaseshift = -np.radians(phase)                           
                return alpha_initial + amp*np.sin(-per*np.arange(num_step) + phaseshift)
            
            case 'custom':
                if np.shape(a1) == (num_step,):
                    return -np.radians(a1)
                else:
                    print("\nImproper custom pitching profile array size\n")
            case _:
                print("\nImproper pitching type selected\n")
    except TypeError:
        print("\nArguments missing or improperly entered in 'pitching' function:\n")

    
alpha_step = pitching('periodic',a1=5,amplitude=15,period=5)           # Call statement for pitching function

"""
Examples of different arguments for 'pitching':
    pitching('periodic',a1=15,amplitude=45,period=5,phase=180)
    pitching('rate',rate=30)
    pitching('constant')
    pitching('custom',np.linspace(0,90,num_step))
"""


# Displacement due to pitching
body_xrot = np.zeros((num_step,num_panels+1))
body_yrot = np.zeros((num_step,num_panels+1))
for i in range(num_step):
    body_xrot[i] = xinitial_body * np.cos(alpha_step[i]) - yinitial_body * np.sin(alpha_step[i])   
    body_yrot[i] = xinitial_body * np.sin(alpha_step[i]) + yinitial_body * np.cos(alpha_step[i])


# Preallocation for time step loop
x_body = np.zeros((num_step,num_panels+1))
y_body = np.zeros((num_step,num_panels+1))
x_vortex = np.zeros((num_step,num_panels))          
y_vortex = np.zeros((num_step,num_panels))
x_colloc = np.zeros((num_step,num_panels))
y_colloc = np.zeros((num_step,num_panels))
normal_body = np.ones((2,num_panels))
motion_body = np.ones((2,num_panels))

x_wake = np.zeros((num_step,num_step))
y_wake = np.zeros((num_step,num_step))
gamma_w = np.zeros((num_step,num_step))
gamma_p = np.zeros((2,num_panels))
circulation = np.zeros((num_step,num_panels+1))
lift = np.ones((num_step))
coeff_lift = np.ones((num_step))

progress = "["+50*"_"+"]     "

for i in range(num_step):
    # Total displacement of the body
    x_body[i,:] = body_xrot[i] + u_disp[i]
    y_body[i,:] = body_yrot[i]

    # Build vortex points
    x_vortex[i,:] = (x_body[i,1:] - x_body[i,0:-1]) * 0.25 + x_body[i,0:-1]
    y_vortex[i,:] = (y_body[i,1:] - y_body[i,0:-1]) * 0.25 + y_body[i,0:-1]

    # Build collocation points
    x_colloc[i,:] = (x_body[i,1:] - x_body[i,0:-1]) * 0.75 + x_body[i,0:-1]
    y_colloc[i,:] = (y_body[i,1:] - y_body[i,0:-1]) * 0.75 + y_body[i,0:-1]

    # Build vectors normal to the body
    normal_body[0,:] = (y_colloc[i] - y_vortex[i])/np.sqrt((y_colloc[i] - y_vortex[i])**2 + (x_colloc[i] - x_vortex[i])**2)
    normal_body[1,:] = -(x_colloc[i] - x_vortex[i])/np.sqrt((y_colloc[i] - y_vortex[i])**2 + (x_colloc[i] - x_vortex[i])**2)

    # Build vectors for body motion
    motion_body[0,:] = u_step[i]
    motion_body[1,:] = 0

    # Update location of trailing edge vorticies
    if i > 0:
        gamma_w[i,1:i+1] = gamma_w[i-1,0:i]
        x_wake[i,1:i+1] = x_wake[i-1,0:i]
        y_wake[i,1:i+1] = y_wake[i-1,0:i]
    x_wake[i,0] = (x_body[i,-1] - x_colloc[i,-1])                                                         \
        /np.sqrt((x_colloc[i,-1] - x_body[i,-1])**2 + (y_colloc[i,-1] - y_body[i,-1])**2)               \
        * delT*abs(u_step[i]) + x_body[i,-1]
    y_wake[i,0] = (y_body[i,-1] - y_colloc[i,-1])                                                         \
        /np.sqrt((x_colloc[i,-1] - x_body[i,-1])**2 + (y_colloc[i,-1] - y_body[i,-1])**2)               \
        * delT*abs(u_step[i]) + y_body[i,-1]

    # Preallocation for loop over panels
    wake_points = i+1
    induced_vel_coeff = np.ones((num_panels,num_panels,2))
    wake_induced_vel_coeff = np.ones((num_panels,wake_points,2))
    
    A = np.ones((num_panels+1,num_panels+1))
    B = np.ones((num_panels+1))

    for j in range(num_panels):

        # Build induced velocity coefficient matrix
        for k in range(num_panels):
            dist_ind_vel = np.sqrt((x_colloc[i,j] - x_vortex[i,k])**2 + (y_colloc[i,j] - y_vortex[i,k])**2)
            induced_vel_coeff[j,k,0] = 1/2/np.pi / dist_ind_vel * (y_colloc[i,j] - y_vortex[i,k]) / dist_ind_vel
            induced_vel_coeff[j,k,1] = 1/2/np.pi / dist_ind_vel * (-x_colloc[i,j] + x_vortex[i,k]) / dist_ind_vel
            """
            induced_vel_coeff[j,k,0] = 1/2/np.pi / np.sqrt((x_colloc[i,j] - x_vortex[i,k])**2 + (y_colloc[i,j] - y_vortex[i,k])**2)     \
                  * (y_colloc[i,j] - y_vortex[i,k]) / np.sqrt((x_colloc[i,j] - x_vortex[i,k])**2 + (y_colloc[i,j] - y_vortex[i,k])**2)
            induced_vel_coeff[j,k,1] = 1/2/np.pi / np.sqrt((x_colloc[i,j] - x_vortex[i,k])**2 + (y_colloc[i,j] - y_vortex[i,k])**2)     \
                  * (-x_colloc[i,j] + x_vortex[i,k]) / np.sqrt((x_colloc[i,j] - x_vortex[i,k])**2 + (y_colloc[i,j] - y_vortex[i,k])**2)
            """
        for k in range(wake_points):
            dist_ind_vel_wake = np.sqrt((x_colloc[i,j] - x_wake[i,k])**2 + (y_colloc[i,j] - y_wake[i,k])**2)
            wake_induced_vel_coeff[j,k,0] = 1/2/np.pi / dist_ind_vel_wake * (y_colloc[i,j] - y_wake[i,k]) / dist_ind_vel_wake
            wake_induced_vel_coeff[j,k,1] = 1/2/np.pi / dist_ind_vel_wake * (-x_colloc[i,j] + x_wake[i,k]) / dist_ind_vel_wake
            """
            wake_induced_vel_coeff[j,k,0] = 1/2/np.pi / np.sqrt((x_colloc[i,j] - x_wake[k])**2 + (y_colloc[i,j] - y_wake[k])**2)        \
                  * (y_colloc[i,j] - y_wake[k]) / np.sqrt((x_colloc[i,j] - x_wake[k])**2 + (y_colloc[i,j] - y_wake[k])**2)
            wake_induced_vel_coeff[j,k,1] = 1/2/np.pi / np.sqrt((x_colloc[i,j] - x_wake[k])**2 + (y_colloc[i,j] - y_wake[k])**2)        \
                  * (-x_colloc[i,j] + x_wake[k]) / np.sqrt((x_colloc[i,j] - x_wake[k])**2 + (y_colloc[i,j] - y_wake[k])**2)
            """
        
        # Build A and B matrices
        A[j,:-1] = induced_vel_coeff[j,:,0] * normal_body[0,j] + induced_vel_coeff[j,:,1] * normal_body[1,j]
        B[j] = motion_body[0,j] * normal_body[0,j] + motion_body[1,j] * normal_body[1,j]                                                \
                - sum(wake_induced_vel_coeff[j,1:,0] * gamma_w[i,1:i+1] * normal_body[0,j]                                                \
                + wake_induced_vel_coeff[j,1:,1] * gamma_w[i,1:i+1] * normal_body[1,j])                                       
        A[j,num_panels] = wake_induced_vel_coeff[j,0,0] * normal_body[0,j] + wake_induced_vel_coeff[j,0,1] * normal_body[1,j]
    
    A[-1,:] = 1
    B[-1] = -sum(gamma_w[i,1:i+1])

    # Solve circulation 
    circulation[i,:] = np.linalg.solve(A,B)

    gamma_p[1,:] = gamma_p[0,:]
    gamma_p[0,:] = circulation[i,0:-1]

    gamma_w[i,0] = circulation[i,-1]

    # Solve Lift and Coeff of Lift
    steadyterm = 0
    unsteadyterm = 0

    for j in range(num_panels):                                             # loop over panels after circulation is calculated
        steadyterm = density*abs(u_step[i])*gamma_p[0,j] + steadyterm
        unsteadyterm = density*chord/num_panels * (sum(gamma_p[0,:j+1]) - sum(gamma_p[1,:j+1]))/delT + unsteadyterm

    lift[i] = steadyterm + unsteadyterm
    coeff_lift[i] = lift[i] / (0.5*density*chord*u_step[i]**2)

    # Wake Roll-up
    rollup_matrix = np.ones((wake_points,wake_points+num_panels,2))
    x_total = np.append(x_vortex[i,:],x_wake[i,:wake_points])
    y_total = np.append(y_vortex[i,:],y_wake[i,:wake_points])
    gamma_total = np.append(np.transpose(circulation[i,:num_panels]),gamma_w[i,:wake_points])


    for j in range(wake_points):                        # for each wake point that currently exists
        for k in range(wake_points+num_panels):             # compare to all other wake points and vortex points
            dist_rollup = np.sqrt((x_wake[i,j] - x_total[k])**2 + (y_wake[i,j] - y_total[k])**2 + (0.04*chord)**2)
            rollup_matrix[j,k,0] = 1/2/np.pi / dist_rollup * (y_wake[i,j] - y_total[k]) / dist_rollup
            rollup_matrix[j,k,1] = 1/2/np.pi / dist_rollup * (-x_wake[i,j] + x_total[k]) / dist_rollup
            
            """
            rollup_matrix[j,k,0] = 1/2/np.pi / np.sqrt((x_wake[j] - x_total[k])**2 + (y_wake[j] - y_total[k])**2 + (0.04*chord)**2)      \
                  * (y_wake[j] - y_total[k]) / np.sqrt((x_wake[j] - x_total[k])**2 + (y_wake[j] - y_total[k])**2 + (0.04*chord)**2)
            rollup_matrix[j,k,1] = 1/2/np.pi / np.sqrt((x_wake[j] - x_total[k])**2 + (y_wake[j] - y_total[k])**2 + (0.04*chord)**2)      \
                  * (-x_wake[j] + x_total[k]) / np.sqrt((x_wake[j] - x_total[k])**2 + (y_wake[j] - y_total[k])**2 + (0.04*chord)**2)
            """

        x_wake[i,j] = x_wake[i,j] + delT*sum(rollup_matrix[j,:,0]*gamma_total)
        y_wake[i,j] = y_wake[i,j] + delT*sum(rollup_matrix[j,:,1]*gamma_total)


    # Progress Bar
    if (i+1) % (num_step//50) == 0:
        print("Progress: ",progress,end="",flush=True)
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
xaxis_size = xaxis_max - xaxis_min

yaxis_min = -1          # Current configuration has y-axis prescribed. Could be made dynamic using same method as x-axis
yaxis_max = 1
yaxis_size = yaxis_max - yaxis_min

axis_ratio = xaxis_size/yaxis_size
plotscalefactor = 3     # Used to scale plot to be a reasonable size. Lower if plot gets too big to prevent graphical warping

# Plot Setup
fig, ax = plt.subplots(figsize=(axis_ratio*plotscalefactor,plotscalefactor))    
    # figsize calculated to ensure consistent scale between axes and prevent graphical distortion of bodies
bodyplot, = ax.plot([], [], 'k', animated=True)
xdata_body, ydata_body = [], []

vortexplot, = ax.plot([], [], 'r.', animated=True)
xdata_vortex, ydata_vortex = [], []

collocplot, = ax.plot([], [], 'c.', animated=True)
xdata_colloc, ydata_colloc = [], []


wakeplot = ax.scatter([], [],s=10,c=[],cmap='RdBu_r',animated=True,vmax=gamma_w.flat[abs(gamma_w).argmax()],vmin=-gamma_w.flat[abs(gamma_w).argmax()])
xdata_wake, ydata_wake = [], []
colorbar = plt.colorbar(wakeplot)


# Toggle visibility of different elements in the animated plot
display_body = True
display_vortex = False
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

    # Vortex Point Plotting
    xdata_vortex = x_vortex[frame]
    ydata_vortex = y_vortex[frame]
    vortexplot.set_data(xdata_vortex, ydata_vortex)

    # Collocation Point Plotting
    xdata_colloc = x_colloc[frame]
    ydata_colloc = y_colloc[frame]
    collocplot.set_data(xdata_colloc,ydata_colloc)

    # Wake Point Plotting
    xdata_wake = np.flip(x_wake,1)[frame,-frame:]
    ydata_wake = np.flip(y_wake,1)[frame,-frame:]
    data = np.stack([xdata_wake,ydata_wake]).T
    wakeplot.set_offsets(data)
    wakeplot.set_array(np.flip(gamma_w,1)[frame,-frame:])
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

end_time = time.process_time()
print("\r"+" "*100,end="",flush=True)
print("\r", end="", flush=True)
print("\nProcess time: ",end_time-start_time," s\n")

# Animation Display
ani = FuncAnimation(fig, update, frames=num_step, init_func=init, blit=True, interval=10, repeat_delay=1000)
    # 10ms between frames, 1000ms pause between playthroughs
plt.show()

