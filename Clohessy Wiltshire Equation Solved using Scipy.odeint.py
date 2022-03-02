#!/usr/bin/env python
# coding: utf-8

# In[14]:


# Import the required modules
import numpy as np
import matplotlib.pyplot as plt
# This makes the plots appear inside the notebook
# get_ipython().run_line_magic('matplotlib', 'inline')


# In[90]:


from scipy.integrate import odeint

miu = 3.986*10**14
a = 6793137
n = np.sqrt(miu/(a**3))
# n = 10
# Define a function which calculates the derivative
def dU_dx(U, x):
    return [U[1],   3*(n**2)*U[0]+2*n*U[3],   U[3],   -2*n*U[1],     U[5],      -(n**2)*U[4]]

# U0 = [x vx y vy z vz]
# NMT IC: [x, 0, 0, -2*n*x, z, vz]
U0 = [1,0,0,-2*n,0,0]
xs = np.linspace(0,10000,1000)

Us = odeint(dU_dx, U0, xs)
xval = Us[:,0] #first column x-coordinates
yval = Us[:,2]
zval = Us[:,4]
print(Us)


# In[91]:


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_aspect('auto')
ax.plot3D(xval,yval,zval, 'grey')

plt.show()

# In[82]:


# Plot the numerical solution
# plt.rcParams.update({'font.size': 14})  # increase the font size
# plt.xlabel("x")
# plt.ylabel("y")
# plt.plot(xval, yval);


# In[70]:


# plt.plot(yval,zval);


# In[77]:



# import numpy as np
# import matplotlib.pyplot as plt
#
# theta = np.linspace( 0 , 2 * np.pi , 150 )
#
# radius = 0.4
#
# a = radius * np.cos( theta )
# b = radius * np.sin( theta )
#
# figure, axes = plt.subplots( 1 )
#
# axes.plot( a, b )
# axes.set_aspect( 1 )
#
# plt.show()
