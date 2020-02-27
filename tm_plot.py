import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

local_t = np.linspace(0 , 5.00000000000001e-2, 100)
local_var_1 = np.linspace(-1 , 1, 100)
local_var_2 = np.linspace(-1 , 1, 100)
local_var_3 = np.linspace(-1 , 1, 100)
local_var_4 = np.linspace(-1 , 1, 100)

xvar =  np.linspace(-2.50626520194835e-6 , 2.50626520200080e-6,  100)
yvar =  np.linspace(-2.50630999981459e-5 , 2.50630999981511e-5, 100)
x3var =  np.linspace(-5.00000000000001e-4 , 5.00000000000001e-4, 100)

#print (local_t)
#print (tvar)


#x = 9.98750260394967e-1 * local_t -2.49895846353391e-3 * local_t * local_var_3 -1.24843782549370e-3 * local_t * local_var_3**2 + 1.04123269313914e-6* local_t *  local_var_3**3 + 2.60091213644523e-7 * local_t * local_var_3**4 -1.30154086642391e-10 * local_t * local_var_3**5  -2.16742678037102e-11 * local_t * local_var_3**6 + 7.74726706204712e-15 * local_t * local_var_3**7 + xvar
#y = 4.99791692706784e-2 * local_t + 4.99375130197484e-2 * local_t * local_var_3  -6.24739615883479e-5 * local_t * local_var_3**2 -2.08072970915618e-5 * local_t * local_var_3**3 + 1.30154086642392e-8 * local_t * local_var_3**4 +2.60091213644523e-9 * local_t * local_var_3**5 -1.08461738868659e-12 * local_t * local_var_3**6 -1.54816198597930e-13 * local_t * local_var_3**7 +    yvar
#x3 = 5.00000000000001e-2 + 5.00000000000000e-2  * local_var_3 + x3var


x = 5.00000000000000e-1 + 5.00000000000000e-1 * local_var_1 + 9.98750260394967e-1 * local_t -2.49895846353391e-3 * local_t * local_var_3 -1.24843782549370e-3 * local_t * local_var_3**2 + 1.04123269313914e-6 * local_t * local_var_3**3 + 2.60091213644523e-7* local_t * local_var_3**4 -1.30154086642391e-10 * local_t * local_var_3**5 -2.16742678037102e-11 * local_t * local_var_3**6 + 7.74726706204712e-15* local_t * local_var_3**7 + xvar

y = 7.50000000000000e-1 + 2.50000000000000e-1* local_var_2 + 4.99791692706784e-2 * local_t +  4.99375130197484e-2 * local_t * local_var_3 -6.24739615883479e-5 * local_t * local_var_3**2  -2.08072970915618e-5 * local_t * local_var_3**3 +1.30154086642392e-8 * local_t * local_var_3**4 +2.60091213644523e-9 * local_t * local_var_3**5 -1.08461738868659e-12 * local_t * local_var_3**6 -1.54816198597930e-13 * local_t * local_var_3**7 + yvar

x3 = 5.00000000000000e-2  +  5.00000000000001e-2 * local_var_3 + x3var


#print (x)
#print (y)
#print (x3)

fig = plt.figure()

ax = fig.gca(projection='3d')
ax.set_xlabel('x1')
ax.set_ylabel('x2')
ax.set_zlabel('x3')
ax.scatter3D(x, y, x3, zdir='z',cmap='viridis')
plt.show()

plt.plot(x, y, 'o', color='black');
plt.show()
