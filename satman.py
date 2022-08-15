#%%
from astropy import units as u
from astropy import time
import numpy as np
from scipy.integrate import odeint
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
from poliastro.util import time_range
from poliastro.plotting import OrbitPlotter3D
from poliastro.plotting import OrbitPlotter2D
import matplotlib.pyplot as plt
import matplotlib.animation as animation
%matplotlib ipympl

# Define Constants
miu = 3.986*10**14
rad_tar = 6378.1+1020 #radius of target from earth surface(km), poliastro uses Earth radius of 6378km
n = np.sqrt(miu/(rad_tar**3))
mass = 1.33 #(Kg) for cubesat 
G = 6.67430*10**-11 #Gravitational constant
mass_E = 5.972*10**24 #Earth Mass
meshsize = 1000

circularorbit = Orbit.circular(Earth, (rad_tar-6378.1) * u.km, 0*u.deg, 0*u.deg, 0*u.deg)#, time.Time("2006-01-19", scale='utc'))inc=0 * u.deg, raan=0 * u.deg, arglat=0 * u.deg,
coord = circularorbit.sample(meshsize)

plotter = OrbitPlotter3D()
plotter.set_attractor(Earth)
plotter.plot(circularorbit)
plotter.set_view(30 * u.deg, 260 * u.deg, distance=3 * u.km)


#Target Satellite Velocity: dependent on distance from Earth
def vel_tar(rad_tar):
    return np.sqrt(G*mass_E/rad_tar)

def controller(): # Move away test
    return np.array([0,0,0]) 
    
def model(U,x): #Add contoller 3/23
    F = controller()
    return [U[3],
            U[4],
            U[5], 
            3*(n**2)*U[0]+2*n*U[4] + F[0]/mass, 
            -2*n*U[3] + F[1]/mass, 
            -(n**2)*U[2] + F[2]/mass]

def run_model(q0):
    q = q0
    
    # length of one circular orbital period for target
    circum = 2*np.pi*rad_tar
    # time steps calculated for 1 orbital period, partitioned in 1000 steps
    t_step = np.linspace(0,circum/vel_tar(rad_tar),meshsize)
    
    Us = odeint(model, q, t_step)
    
    return Us[:,0],Us[:,1],Us[:,2]
    # xval = Us[:,0] #first column x-coordinates


# Initial Condition
U0 = [500,0,0, 0,-2*n*500,0]

rel_x,rel_y,rel_z = run_model(U0)

chaser_x = rel_x + np.copy(coord.x.value)
chaser_y = rel_y + np.copy(coord.y.value)
chaser_z = rel_z + np.copy(coord.z.value)
  
chaser_x = chaser_x*u.km    
chaser_y = chaser_y*u.km 
chaser_z = chaser_z*u.km

angles = np.linspace(0,2*np.pi,1000)

############ Relative Distance Plots  ######################################################
#plot 1:
plt.subplot(1, 3, 1)
plt.plot(angles,rel_x) 
plt.xlabel("") 
plt.ylabel("Relative x distance") 

#plot 2:
plt.subplot(1, 3, 2)
plt.plot(angles,rel_y) 
plt.xlabel("") 
plt.ylabel("Relative y distance")

#plot 3:
plt.subplot(1, 3, 3)
plt.plot(angles,rel_z) 
plt.xlabel("") 
plt.ylabel("Relative z distance")



############ Animated Trajectory ######################################################

#Create Earth on Figure 
u = np.linspace(0,2*np.pi,100)
v = np.linspace(0,  np.pi,100)
r = 6378.1 #km
xe = r * np.outer(np.cos(u), np.sin(v))
ye = r * np.outer(np.sin(u), np.sin(v))
ze = r * np.outer(np.ones(np.size(u)), np.cos(v))

#PRESET CIRCULAR POINTS for target ##############################
x = coord.x.value
y = coord.y.value
z = coord.z.value

#Put into array of shape(number of values, 3)
coordinate = x.reshape(x.shape[0],1)
y = y.reshape(x.shape[0],1)
z = z.reshape(x.shape[0],1)
coordinate= np.append(coordinate, y, axis = 1)
coordinate= np.append(coordinate, z, axis = 1)
#PRESET POINTS for chaser ##############################
xc = chaser_x.value
yc = chaser_y.value
zc = chaser_z.value

#Put into array of shape(number of values, 3)
coordinatec = xc.reshape(xc.shape[0],1)
yc = yc.reshape(xc.shape[0],1)
coordinatec= np.append(coordinatec, yc, axis = 1)
zc = zc.reshape(xc.shape[0],1)
coordinatec= np.append(coordinatec, zc, axis = 1)

# Animation ############################################

def update_lines(num, walks, lines):
    for line, walk in zip(lines, walks):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(walk[:num, :2].T)
        line.set_3d_properties(walk[:num, 2])
    return lines

num_steps = x.shape[0]
walks = [coordinate]
walksc = [coordinatec]


# Attaching 3D axis to the figure
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.plot_surface(xe,ye,ze,rstride = 5, cstride = 5)

# Create lines initially without data
lines = [ax.plot([], [], [])[0]]
linesc = [ax.plot([], [], [])[0]]

# Setting the axes properties
ax.set(xlim3d=(-15000, 15000), xlabel='X')
ax.set(ylim3d=(-15000, 15000), ylabel='Y')
ax.set(zlim3d=(-15000, 15000), zlabel='Z')

# Creating the Animation object
ani = animation.FuncAnimation(
    fig, update_lines, num_steps, fargs=(walks, lines), interval=2)
anic = animation.FuncAnimation(
    fig, update_lines, num_steps, fargs=(walksc, linesc), interval=2)
plt.show()

############ Wire Frame Plot Trajectory ######################################################
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z') 
u = np.linspace(0,2*np.pi,20)
v = np.linspace(0,  np.pi,20)
r = 6378.1 #km
xe = r * np.outer(np.cos(u), np.sin(v))
ye = r * np.outer(np.sin(u), np.sin(v))
ze = r * np.outer(np.ones(np.size(u)), np.cos(v)) 
ax.set_aspect('auto')
ax.set(xlim3d=(-15000, 15000), xlabel='X')
ax.set(ylim3d=(-15000, 15000), ylabel='Y')
ax.set(zlim3d=(-15000, 15000), zlabel='Z')
ax.plot_wireframe(xe,ye,ze,color='grey',linewidth = 0.4)
ax.plot3D(coord.x.value,coord.y.value,coord.z.value, 'green',linewidth = 1.2, label = "Target Satellite Orbit")
ax.plot3D(chaser_x,chaser_y,chaser_z,color = 'red',linewidth = 1.2, label = "Chaser Satellite Orbit")
ax.legend(loc='upper left')
ax.view_init(60, 35)
plt.show()
# %%


