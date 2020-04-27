import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

local_t = np.linspace(0 , 5.00000000000001e-2, 100)
local_var_1 = np.linspace(-1 , 1, 100)
local_var_2 = np.linspace(-1 , 1, 100)
local_var_3 = np.linspace(-1 , 1, 100)
local_var_4 = np.linspace(-1 , 1, 100)

#xvar =  np.linspace(-2.50626520194835e-6 , 2.50626520200080e-6,  100)
#yvar =  np.linspace(-2.50630999981459e-5 , 2.50630999981511e-5, 100)
#x3var =  np.linspace(-5.00000000000001e-4 , 5.00000000000001e-4, 100)

#print (local_t)
#print (tvar)


#x = 9.98750260394967e-1 * local_t -2.49895846353391e-3 * local_t * local_var_3 -1.24843782549370e-3 * local_t * local_var_3**2 + 1.04123269313914e-6* local_t *  local_var_3**3 + 2.60091213644523e-7 * local_t * local_var_3**4 -1.30154086642391e-10 * local_t * local_var_3**5  -2.16742678037102e-11 * local_t * local_var_3**6 + 7.74726706204712e-15 * local_t * local_var_3**7 + xvar
#y = 4.99791692706784e-2 * local_t + 4.99375130197484e-2 * local_t * local_var_3  -6.24739615883479e-5 * local_t * local_var_3**2 -2.08072970915618e-5 * local_t * local_var_3**3 + 1.30154086642392e-8 * local_t * local_var_3**4 +2.60091213644523e-9 * local_t * local_var_3**5 -1.08461738868659e-12 * local_t * local_var_3**6 -1.54816198597930e-13 * local_t * local_var_3**7 +    yvar
#x3 = 5.00000000000001e-2 + 5.00000000000000e-2  * local_var_3 + x3var


#x = [5.00000000000000e-1 , 5.00000000000000e-1] + 5.00000000000000e-1 * local_var_1 + 9.98750260394967e-1 * local_t -2.49895846353391e-3 * local_t * local_var_3 -1.24843782549370e-3 * local_t * local_var_3**2 + 1.04123269313914e-6 * local_t * local_var_3**3 + 2.60091213644523e-7* local_t * local_var_3**4 -1.30154086642391e-10 * local_t * local_var_3**5 -2.16742678037102e-11 * local_t * local_var_3**6 + 7.74726706204712e-15* local_t * local_var_3**7 + xvar
#y = 7.50000000000000e-1 + 2.50000000000000e-1* local_var_2 + 4.99791692706784e-2 * local_t +  4.99375130197484e-2 * local_t * local_var_3 -6.24739615883479e-5 * local_t * local_var_3**2  -2.08072970915618e-5 * local_t * local_var_3**3 +1.30154086642392e-8 * local_t * local_var_3**4 +2.60091213644523e-9 * local_t * local_var_3**5 -1.08461738868659e-12 * local_t * local_var_3**6 -1.54816198597930e-13 * local_t * local_var_3**7 + yvar
#x3 = 5.00000000000000e-2  +  5.00000000000001e-2 * local_var_3 + x3var

x = min(2.50000000000000e-1 , 2.50000000000000e-1) + min(2.50000000000000e-1 , 2.50000000000000e-1) * local_var_1 + min(4.79830826328700 , 4.79830826328701) * local_t + min(-4.00649368095794e-1 , -4.00649368095793e-1) * local_t * local_var_3 + min(-1.94871294342744e-1 , -1.94871294342743e-1) * local_t * local_var_3**2 + min(5.42379082059680e-3 , 5.42379082059681e-3) * local_t * local_var_3**3 + min(1.31903507358244e-3 , 1.31903507358245e-3) * local_t * local_var_3**4 + min(-2.20273704701488e-5 , -2.20273704701487e-5) * local_t * local_var_3**5 + min(-3.57128746172447e-6 , -3.57128746172446e-6) * local_t * local_var_3**6 + min(4.25993611056627e-8 , 4.25993611056628e-8) * local_t * local_var_3**7 + min(5.17996114426017e-9 , 5.17996114426018e-9) * local_t * local_var_3**8 + min(-4.80574042473258e-11 , -4.80574042473257e-11) * local_t * local_var_3**9 + min(-4.67491493269481e-12 , -4.67491493269480e-12) * local_t * local_var_3**10 + min(3.54860241817184e-14 , 3.54860241817185e-14) * local_t * local_var_3**11 + min(-2.85207401672995e-4 , 2.85207401673283e-4)

y = min(3.10000000000000 , 3.10000000000001) + min(3.99999999999999e-1 , 4.00000000000000e-1) * local_var_2 + min(1.40578725647646 , 1.40578725647647) * local_t + min(1.36751785503679 , 1.36751785503680) * local_t * local_var_3 + min(-5.70925349536507e-2 , -5.70925349536506e-2) * local_t * local_var_3**2 + min(-1.85127729625607e-2 , -1.85127729625606e-2) * local_t * local_var_3**3 + min(3.86445095967522e-4 , 3.86445095967523e-4) * local_t * local_var_3**4 + min(7.51849991941993e-5 , 7.51849991941994e-5) * local_t * local_var_3**5 + min(-1.04630009733207e-6 , -1.04630009733206e-6) * local_t * local_var_3**6 + min(-1.45402418084497e-7 , -1.45402418084496e-7) * local_t * local_var_3**7 + min(1.51760223938923e-9 , 1.51760223938924e-9) * local_t * local_var_3**8 + min(1.64032102901572e-10 , 1.64032102901573e-10) * local_t * local_var_3**9 + min(-1.36963602104879e-12 , -1.36963602104878e-12) * local_t * local_var_3**10 + min(-1.21122795983457e-13 , -1.21122795983456e-13) * local_t * local_var_3**11 + min(-5.40202298494961e-4 , 5.40202298495045e-4)

x3 = min(2.84999999999999e-1 , 2.85000000000000e-1) + min(2.84999999999999e-1 , 2.85000000000000e-1) * local_var_3 + min(-1.00000000000001e-3 , 1.00000000000001e-3)

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
