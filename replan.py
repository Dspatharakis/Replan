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
    condition1 = True
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

    print (x1in,x2in,x3in,d,xd1,xd2,xd3,condition1,condition2)


if __name__ == '__main__':
    simulation()

