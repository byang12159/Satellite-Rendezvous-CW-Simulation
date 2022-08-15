from astropy import units as astro
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

class CWHmotion:
    def __init__(self, miu = 3.986*10**14, rad_tar = 6378.1+1020 #radius of target from earth surface(km), poliastro uses Earth radius of 6378km,
                 ,mass = 1.33 #(Kg) for cubesat 
                 ,G = 6.67430*10**-11 #Gravitational constant
                 ,mass_E = 5.972*10**24 #Earth Mass
                 ,meshsize = 1000, periods = 1
                 ,initialstate = np.zeros(6)
                 ,thrustforce = np.array([0,0,0])
                 ):
        
        ''' System Parameters'''
        self.miu = miu
        self.rad_tar = rad_tar
        self.mass = mass
        self.G = G
        self.mass_E = mass_E
        self.meshsize = meshsize
        self.peiods = periods
        self.thrustforce = thrustforce
        
        self.vel_tar = np.sqrt(G*mass_E/rad_tar)
        self.n = np.sqrt(miu/(rad_tar**3)) #orbital rate of target body, see wiki
        self.initialstate = initialstate

    
    def main(self):
        
        rel_x,rel_y,rel_z = self.runmodel()
        #radius of target from earth surface(km), poliastro uses Earth radius of 6378km,
        circularorbit = Orbit.circular(Earth, (sim.rad_tar-6378.1) * astro.km, 0*astro.deg, 0*astro.deg, 0*astro.deg)
        coord = circularorbit.sample(sim.meshsize)
        
        chaser_x = rel_x + np.copy(coord.x.value)
        chaser_y = rel_y + np.copy(coord.y.value)
        chaser_z = rel_z + np.copy(coord.z.value)

        # # ################################################################ Relative Distance Plots  ######################################################
        
        angles = np.linspace(0,2*np.pi,1000)
        
        return rel_x,rel_y,rel_z, angles, chaser_x, chaser_y, chaser_z, coord
        
    def controller(self):
        return [self.thrustforce[0]/self.mass, self.thrustforce[1]/self.mass, self.thrustforce[2]/self.mass] 
    
    def reltrajectorygenerate(self,U,x):
        accel = self.controller()
        return [U[3], 
                U[4], 
                U[5], 
                3*(self.n**2)*U[0]+2*self.n*U[4] + accel[0], 
                -2*self.n*U[3] + accel[1], 
                -(self.n**2)*U[2] + accel[2]]
    
    def runmodel(self):
        q = self.initialstate
    
        # length of one circular orbital period for target
        circum = 2*np.pi*self.rad_tar
        
        # time steps calculated for 1 orbital period, partitioned in 1000 steps
        t_step = np.linspace(0,circum/self.vel_tar,self.meshsize)

        Us = odeint(self.reltrajectorygenerate, q, t_step)

        return Us[:,0],Us[:,1],Us[:,2]
        # xval = Us[:,0] #first column x-coordinates


if __name__ == '__main__':

    # Take input as result from multi_sat.py, the last entries of X_NMT 
    # X_NMT_last = [X_NMT[0,-1],X_NMT[1,-1],X_NMT[2,-1],X_NMT[3,-1],X_NMT[4,-1],X_NMT[5,-1]]
        
    X_NMT_last = np.array([-1000.1369261894047, -1484.113395146748, 0.0, -0.7620922284079551, 2.054281246393022, 0.0])
    sim = CWHmotion()
    sim.initialstate = X_NMT_last
    # sim.thrustforce = np.array([1,2,3])
    
    rel_x,rel_y,rel_z, angles, chaser_x, chaser_y, chaser_z, coord = sim.main()
    
    
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
    plt.show()


    # # ############ Wire Frame Plot Trajectory ######################################################
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
    # ax.plot3D(X_NMT[0,:],X_NMT[1,:],X_NMT[2,:])
    # ax.plot3D(Y_NMT[0,:],Y_NMT[1,:],Y_NMT[2,:])
    ax.legend(loc='upper left')
    ax.view_init(60, 35)
    plt.show()
    
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
    xc = chaser_x
    yc = chaser_y
    zc = chaser_z

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

    print("complete")
 