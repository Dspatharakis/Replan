#!/usr/bin/env python3
import cvxpy as cp
import numpy as np
import math

# v3 = cp.Variable()
# v2 = cp.Variable()
# z3 = cp.Variable()
# z2 = cp.Variable()

# # DCP problems.
# prob1 = cp.Problem(cp.Minimize(cp.square(v3-z3)),
#                     [cp.square(v3 - z3)   <= 5  , 
#                     #v3 - z3   <= 2.5 , 
#                     #v1 - z1   <= -cp.sqrt(2.5 - cp.square(v2 - z2)) , 
#                     z3 >= 0 , z3 <= 1, 
#                     v3 >= 2 , v3 <= 4, ])


# try:
#     print(prob1.solve())
# except Exception as e:
#     print(e)

# prob1 = cp.Problem(cp.Maximize(v3),
#                     [v3 - z3   >= 2.2 , v3 - z3   <= 2.5 , 
#                     #v1 - z1   <= -cp.sqrt(2.5 - cp.square(v2 - z2)) , 
#                     z3 >= 0 , z3 <= 1, 
#                     v3 >= 2 , v3 <= 4, ])


# try:
#     print(prob1.solve())
# except Exception as e:
#     print(e)


# v1 = cp.Variable()
# v2 = cp.Variable()
# z1 = cp.Variable()
# z2 = cp.Variable()

# prob1 = cp.Problem(cp.Minimize(cp.norm(v1-z1)),
#                     [#v1 - z1   <= cp.sqrt(2.5 - cp.square(v2 - z2))  , 
#                     z1 >= 0 , z2 >= 0 , z1 <= 4 , z2 <= 5, 
#                     v1 <= 10 , v2 <= 12, v1 >= 0 , v2 >= 0 ])


# try:
#     print(prob1.solve())
# except Exception as e:
#     print(e)

#Gains:  0.2588810084625552 0.2588810084625552 0
#Flow: 23.17672222020196 23.19083092298054 22.81545083957362 22.82966532302303 -0.9502109953154041 -0.947985120333111
#22.81545083957362 22.82966532302303 23.02550000000075 23.040500000001227 0.06701937654258999 2
#23.17672222020196 23.19083092298054 23.02550000000075 23.040500000001227 0.06701937654258999 2
#Optimization 23.17672222020196 23.19083092298054 23.02550000000075 22.82966532302303 -0.9492109953356229 -0.9489851203128921

#vmin,vmax,zmin,zmax,y,w
vmin = 23.17672222020196
vmax = 23.19083092298054
zmin = 23.02550000000075
zmax = 23.040500000001227
y = 2#0.06701937654258999 
w = 0.1 * y


vmin = 14.87665278997093
vmax =  14.877112448508544
zmin =  15 
zmax = 15 
y = 0.026680500000000027
w=  2


v = cp.Variable()
z = cp.Variable()
prob1 = cp.Problem(cp.Minimize(v ), #cp.square(v-z)),
                    [cp.square(v-z)   <= y+w , 
                    z >= zmin , z <= zmax, 
                    v >= vmin , v <= vmax ])

try:
    minimum = prob1.solve()
except Exception as e:
    print(e)

v = cp.Variable()
z = cp.Variable()
prob1 = cp.Problem(cp.Maximize(v), # v-z),
                    [cp.square(v-z)   <= y-w , 
                    z >= zmin , z <= zmax, 
                    v >= vmin , v <= vmax ])

try:
    maximum = prob1.solve()
except Exception as e:
    print(e)


print (minimum)
print (maximum)