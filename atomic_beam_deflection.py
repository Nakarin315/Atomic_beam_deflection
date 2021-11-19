"""
Created on Thu Oct  21 18:20:54 2021

@author: Nakarin
"""
##############################
#Declarations and Setup
from math import pi,sqrt,floor,fabs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import*
from scipy import constants
from scipy import optimize as opt

#Universal constants
k_B = constants.k # Boltzmann constant
mu_0 = constants.mu_0 #permeability of free space
hbar = constants.hbar #reduced Planck's constant
mu_B = (constants.physical_constants['Bohr magneton'])[0] # Bohr magneton in SI units
a_m = (constants.physical_constants['atomic mass constant'])[0] # Atomic mass constant in kg
c = constants.c #speed of light in vacuum

#Global experimental parameters
T = 450+273.15 # temperature of oven, assuming same parameters as Schioppo (2012)
beam_radius = 0.5 # radius of atomic beam in cm
m_Sr = 86.9088774970*a_m; # atomic mass in kg
w_c = 460.8622e-9 # cooling wavelength in m
k_wave=2*pi/w_c; # wave number
Gamma = 2*pi*32e6 # decay rate of cooling transition in angular frequency (includes factor of 2*pi)
delta = -58e6; # Detuning
dt = 5e-8 # Time step
#Saturation intensity of cooling transition
Isat = 2 * hbar * (pi**2) * c * Gamma / (3 * (w_c**3)) # in SI
I_laser = 0.1*Isat;
s0 = I_laser/Isat;
B = 0.5 ; # Testla/cm (50 G/cm)
w_beam = 10e-3 # Beam waist of laser


# Size of laser
xmin_laser = -30e-3; xmax_laser = 30e-3;
ymin_laser = -w_beam/2; ymax_laser = w_beam/2;
# Boundary of plot
xmin_bound = -50e-3; xmax_bound = 50e-3;
ymin_bound = -80e-3; ymax_bound = 50e-3;

def acceleration_x(x,y,vx,vy): # Acceleration due to MOT beam
        # Gaussian function represents a profile of laser
    return -np.exp(-y**2/(2*w_beam**2))/(2*pi*w_beam)*(hbar*k_wave*s0*Gamma/(2*m_Sr))*((1/(1+s0+4*(delta-k_wave*vx-mu_B*B*x)**2/Gamma**2))-(1/(1+s0+4*(delta+k_wave*vx+mu_B*B*x)**2/Gamma**2)))

def acceleration_y(x,y): # No acceleration along y-axis
    return 0

fig = plt.figure(2)
ax = fig.add_subplot(111)
theta =pi/3; # Angle between atomic beam and laser beam
# Initial position of atomic beam
x_traj = [xmin_bound*np.cos(theta)*8/9,] 
y_traj = [ymin_bound*np.sin(theta)*8/9,]
x_F=np.linspace(xmin_bound,xmax_bound,400);
y_F = np.linspace(ymin_bound,ymax_bound,400);
# Plot laser beam
X_F, Y_F = np.meshgrid(x_F, y_F)
laser_beam = np.exp(-(Y_F**2)/(2*w_beam**2))/(2*pi*w_beam)
plt.contour(X_F,Y_F,laser_beam,2000,cmap=cm.Blues)
plt.xlim([xmin_bound, xmax_bound])
plt.ylim([ymin_bound, ymax_bound])

# Calculation of atomic beam trajectory
for i in range(10):
    v_atom=100*i/10+30; # vary velocity from 30 to 120 m/s
    # Trajectory with laser beam
    x_traj = [xmin_bound*np.cos(theta)*8/9,]
    y_traj = [ymin_bound*np.sin(theta)*8/9,]
    # Trajectory without laser beam
    x_traj_no = [x_traj[0],]
    y_traj_no = [y_traj[0],] 
    theta =pi/3
    x_traj = [xmin_bound*np.cos(theta)*8/9,]
    y_traj = [ymin_bound*np.sin(theta)*8/9,]
    vx_traj = [v_atom*np.cos(theta),]
    vy_traj = [v_atom*np.sin(theta),]
    while x_traj[-1] > xmin_bound and x_traj[-1] < xmax_bound and y_traj[-1] > ymin_bound and y_traj[-1] < ymax_bound :
        vx_update = vx_traj[-1] - dt * acceleration_x(x_traj[-1],y_traj[-1],vx_traj[-1],vy_traj[-1])
        vy_update = vy_traj[-1] - dt * acceleration_y(x_traj[-1],y_traj[-1])
        x_update = x_traj[-1] + dt * vx_traj[-1]
        y_update = y_traj[-1] + dt * vy_traj[-1]
        x_traj += [x_update,]
        y_traj += [y_update,]
        x_traj_no += [x_traj_no[-1]+dt*vx_traj[0],]
        y_traj_no += [y_traj_no[-1]+dt*vy_traj[0],]
        vx_traj += [vx_update,]
        vy_traj += [vy_update,]
    plt.plot(x_traj,y_traj, label='%s m/s' % v_atom)
    plt.plot(x_traj_no,y_traj_no,'r--', alpha=0.1)
plt.legend()
plt.show()




# visualization of laser force
@np.vectorize
def acceleration_x(x,y,vx,vy):
    return -np.exp(-y**2/(2*w_beam**2))/(2*pi*w_beam)*(hbar*k_wave*s0*Gamma/(2*m_Sr))*((1/(1+s0+4*(delta-k_wave*vx-mu_B*B*x)**2/Gamma**2))-(1/(1+s0+4*(delta+k_wave*vx+mu_B*B*x)**2/Gamma**2)))

def acceleration_y(x,y):
    return np.exp(-(x**2+y**2))
F = m_Sr*acceleration_x(X_F,Y_F,0.1,0.1);
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X_F, Y_F, F)
