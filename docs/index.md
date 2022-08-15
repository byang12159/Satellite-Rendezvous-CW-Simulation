## Introduction
This repository provides a satellite rendezvous simulation following the Clohessy Whiltshire Equations (CWH).

### Getting Started

#### Installations
1. Astropy ``` pip install astropy ```
2. NumPy ```pip3 install numpy```
3. SciPy ```pip3 install scipy```
4. Matplotlib ```pip3 install matplotlib```

#### Basic Mechanics
$U_{target}$ = state of target satellite relative to Earth <br>
$U_{chaser}$ = state of chaser satellite relative to Earth

The CWH equations describe the relative motion of a chaser spacecraft in circular/elliptical orbit around a target spacecraft in circular orbit.

To update the state of the chaser at any time point, and obtain results as Earth frame instead of relative frame, requires the vector addition as shown
$$ U_{chaser}(t) = U_{target}(t) + U_{chaser}^{update}(t) $$

INSERT diagram 

The update $U_{chaser}^{update}(t)$ is obtained by solving the following equations using scipy's odeint solver
$$ \dot U = AU + BF $$
Matrix $A$ contains the dynamics given by CWH. Matrix $B$ maps the forces in the correct locations. More details can be found in the mechanics section. 

### Usage
#### System Parameters
Most of the defined variables do not need to be changed for Earth-centric cubesat simulations. The following are variables that will differ from mission to mission:

``` python
rad_tar = 6378.1+1020                          # radius of target from earth center (km)
mass = 1.33                                    # Chaser satellite mass (kg) 
meshsize = 1000                                # Granularity of time step
meshsize = 1000                                # Granularity of time step
periods = 1
initialstate = np.array([x, y, z, Vx, Vy, Vz]) # Initial state of chaser at the start of simulation
thrustforce = np.array([Fx, Fy, Fz])           # Chaser thrust force 
```

#### Example Case 
A satellite (chaser) initiated an orbital transfer to approach the target satellite orbit. The simulation shows if it can maintain within a stable orbit, governed by the CWH equation.
``` python
if __name__ == '__main__':

    # Take input as result from multi_sat.py, the last entries of X_NMT 
    # X_NMT_last = [X_NMT[0,-1],X_NMT[1,-1],X_NMT[2,-1],X_NMT[3,-1],X_NMT[4,-1],X_NMT[5,-1]]
        
    X_NMT_last = np.array([-1000.1369261894047, -1484.113395146748, 0.0, -0.7620922284079551, 2.054281246393022, 0.0])
    sim = CWHmotion()
    sim.initialstate = X_NMT_last

    rel_x,rel_y,rel_z, angles, chaser_x, chaser_y, chaser_z, coord = sim.main()
    
    
```


+ https://en.wikipedia.org/wiki/Clohessy%E2%80%93Wiltshire_equations
