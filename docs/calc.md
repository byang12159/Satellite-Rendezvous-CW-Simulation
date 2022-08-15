#### Mechanics of Calculations
##### Variables
$r$ = radius of target spacecraft's circular orbit <br>
$\mu$ = standard gravitational parameter <br>
$n$ = mean motion of the target spacecraft (1)
$$
n = \sqrt{\frac{ {\mu} }{r^3}}
$$ 


##### Clohesy Whiltshire Equation
$$
\ddot x - 2n\dot y - 3n^2x =  a_{x}
$$

$$
\ddot y + 2n\dot x =  a_{y}
$$

$$
\ddot z + n^2 z =  a_{z}
$$

$$
Thrust F_{x} = m_{Satellite}/ a_{x}
$$

$$
Thrust_{y} = F_{y}/m_{Satellite}
$$

$$
Thrust_{x} = F_{z}/m_{Satellite}
$$

The Clohessy Whiltshire update is rewritten in matrix format as:




$\left[\begin{array}{cc} \dot x \\ \dot y \\  \dot z \\ \dot v_{x} \\  \dot v_{y} \\ \dot v_{z} \end{array}\right] =$
$\left[\begin{array}{cc} 0 & 1 & 0 & 0 & 0 & 0 \\  
                         0 & 0 & 0 & 1 & 0 & 0 \\
                         0 & 0 & 0 & 0 & 0 & 1 \\ 
                         3n^2 & 0 & 0 & 2n & 0 & 0 \\
                         0 & -2n & 0 & 0 & 0 & 0 \\
                         0 & 0 & 0 & 0 & -n^2 & 0 
\end{array}\right] *$
$\left[\begin{array}{cc} x \\  y \\  z \\ v_{x} \\ v_{y} \\ v_{z} \end{array}\right] +$
$\left[\begin{array}{cc} 0 & 0 & 0 \\
                         0 & 0 & 0 \\ 
                         0 & 0 & 0 \\
                         1 & 0 & 0 \\  
                         0 & 1 & 0 \\
                         0 & 0 & 1
\end{array}\right] *$
$\left[\begin{array}{cc} F_{x} \\ F_{y} \\ F_{z} \end{array}\right]$

This is mapped out in the Relative Trajectory function:
``` python
    def reltrajectorygenerate(self,U,x):
        accel = self.controller()
        return [U[3], 
                U[4], 
                U[5], 
                3*(self.n**2)*U[0]+2*self.n*U[4] + accel[0], 
                -2*self.n*U[3] + accel[1], 
                -(self.n**2)*U[2] + accel[2]]
```
