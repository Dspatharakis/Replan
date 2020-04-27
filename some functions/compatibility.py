#!/usr/bin/env python3
import cvxpy as cp
import numpy as np
import math

v3 = cp.Variable()
v2 = cp.Variable()
z3 = cp.Variable()
z2 = cp.Variable()

# DCP problems.
prob1 = cp.Problem(cp.Minimize(cp.square(v3-z3)),
                    [cp.square(v3 - z3)   <= 5  , 
                    #v3 - z3   <= 2.5 , 
                    #v1 - z1   <= -cp.sqrt(2.5 - cp.square(v2 - z2)) , 
                    z3 >= 0 , z3 <= 1, 
                    v3 >= 2 , v3 <= 4, ])


try:
    print(prob1.solve())
except Exception as e:
    print(e)

prob1 = cp.Problem(cp.Maximize(v3),
                    [v3 - z3   >= 2.2 , v3 - z3   <= 2.5 , 
                    #v1 - z1   <= -cp.sqrt(2.5 - cp.square(v2 - z2)) , 
                    z3 >= 0 , z3 <= 1, 
                    v3 >= 2 , v3 <= 4, ])


try:
    print(prob1.solve())
except Exception as e:
    print(e)


v1 = cp.Variable()
v2 = cp.Variable()
z1 = cp.Variable()
z2 = cp.Variable()

prob1 = cp.Problem(cp.Minimize(cp.norm(v1-z1)),
                    [#v1 - z1   <= cp.sqrt(2.5 - cp.square(v2 - z2))  , 
                    z1 >= 0 , z2 >= 0 , z1 <= 4 , z2 <= 5, 
                    v1 <= 10 , v2 <= 12, v1 >= 0 , v2 >= 0 ])


try:
    print(prob1.solve())
except Exception as e:
    print(e)


