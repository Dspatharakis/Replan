#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import random
import math
import time 
from numpy import exp,arange , arctan , sqrt
from scipy import signal
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
from matplotlib import cm
from matplotlib.ticker import MaxNLocator


def euler_forwarding(x1,x2,x3,d3,u):
    x3t = x3 + d3
    x1 = x1 + u * 0.5 * math.sin(x3)
    x2 = x2 + u * 0.5 * math.cos(x3)
    return x1,x2,x3t

def discretization(x1,x2,x3,d3,u):
    x3t = x3 + d3*0.5
    x1 = x1 + u /d3 * (math.sin(x3t) - math.sin(x3))
    x2 = x2 + u /d3 * (math.cos(x3) - math.cos(x3t))
    return x1,x2,x3t 

def build_model(time, d1, d2, d3, u1, u2, u3, invtran, jumptime , invrot, initial_set):
#everything must be strings 

    with open('replan.model', 'w') as the_file:
        the_file.write('hybrid reachability\n''{\n''  state var x,y,x3,t \n'' setting\n''{\n'
        'fixed steps 0.05\n''time '+str(time)+'\n' # time  
        'remainder estimation 1e-2\n''identity precondition\n''gnuplot octagon  x,y\n'
        'fixed orders 8\n''cutoff 1e-15\n''precision 2000\n''output ltv1_test\n''max jumps 1\n''print on\n''}\n')
        the_file.write('modes\n''{\n''tran\n''{ \n''nonpoly ode\n''{\n'"t' = 1\n")
        the_file.write("x3' ="+ str(d1) + "\n" # to d 1
        "x' = "+str(u1)+"*cos(x3) \n"   # to u 1
        "y' = "+str(u2)+"*sin(x3) \n""}\n") # to u 2
        the_file.write("inv\n""{\n"
        ""+invtran[0]+"\n"    # add compatibillity constraints 
        "}\n""}\n""rot\n""{\n""nonpoly ode\n""{\n"
        "x3' = "+str(u1)+"\n"  # to u3
        "x' = [-0.1 , 0.1]\n""y' = [-0.1 , 0.1]\n" # to d1
        "t' = 1\n""}\n""inv \n""{\n" # to d2
        ""+invrot[0]+"\n"  # add compatibillity constraints 
        ""+invrot[1]+" \n" # add compatibillity constraints 
        "}\n""}\n""}\n""jumps\n""{\n""tran -> rot\n"
        "guard { t = "+str(jumptime)+"  }\n" # edw einai h fasoula gia to jump
        "reset { }\n""parallelotope aggregation { }\n""rot -> tran \n""guard {x = 5}\n""reset{}#x' := x - 4.9 }\n""parallelotope aggregation { }\n""}\n""init\n""{\n""tran \n""{\n"
        "x * y <= 2""\n" # in "+initial_set[0]+"\n"  # add initial set constraints 
        "y in "+initial_set[1]+"\n"  # add initial set constraints 
        "x3 in "+initial_set[2]+"\n"   # add initial set constraints  
        "}\n""}\n""}\n"
    )
def point_in_poly(reachability_points,x,y):
    #current convex hull
    points = np.vstack(reachability_points)
    hull = ConvexHull(points)#,qhull_options='QJ')
    cx = np.mean(hull.points[hull.vertices,0])
    cy = np.mean(hull.points[hull.vertices,1])
    xvert1 = hull.points[hull.vertices,0]
    #attemp icnluding point
    temp = reachability_points
    temp.append([x,y])
    points = np.vstack(temp)
    hull = ConvexHull(points)#,qhull_options='QJ')
    cx = np.mean(hull.points[hull.vertices,0])
    xvert2 = hull.points[hull.vertices,0]
    if len(xvert1) != len(xvert2): return False   # false stands for point outside of polytope
    else: return True #true stands for point in polytope 


def simulation():
    time = 0.05 
    d1 = [-0.1,0.1]
    d2 = [-0.1,0.1]
    d3 = [-0.01,0.01]
    u1 = 1
    u2 = 1 
    u3 = 1  
    invtran = ["t<=1"]
    jumptime = 1 #choose 1 for only tran or 0.05 for rot 
    invrot = ["x>=0","y>=0"]
    initial_set = ["[0.0 , 1.0]", "[10.5 , 12.0]", "[0, 0.1]"]
    build_model(time, d1, d2, d3, u1, u2, u3, invtran, jumptime , invrot, initial_set) 
    condition1 = False # True
    condition2 = True

    
    while (condition1 and condition2):
        x1in = x1 = random.uniform(0,2)
        x2in = x2 = random.uniform(0,2) 
        x3in = x3 = random.uniform(0,1.8)
        reachability_points = [[2.856305,2.250760],[2.218212,2.888853],[0.034515,2.888853],[-0.234695,2.619642],[-0.234695,0.217331],[0.088161,-0.105526],[2.587095,-0.105526],[2.856305,0.163685]] 

        d = random.uniform(-0.2,0.2)
        u= 1 
        #xe1,xe2,xe3 = euler_forwarding(x1in,x2in,x3in)
        xd1,xd2,xd3 = discretization(x1in,x2in,x3in,d,u)
        xe1,xe2,xe3 = euler_forwarding(x1in,x2in,x3in,d,u)

        condition1 = point_in_poly(reachability_points,xd1,xd2)
        #print (condition1)
        condition2 = point_in_poly(reachability_points,xd1,xd2)
        #print (condition2)
        #check with flow reachability set 

    #print (x1in,x2in,x3in,d,xd1,xd2,xd3,condition1,condition2)


if __name__ == '__main__':
    simulation()

