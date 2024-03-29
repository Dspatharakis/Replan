#!/usr/bin/env python3
from re import T
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import time 
import csv
import os
import subprocess
import itertools
import cvxpy as cp
from pypoman import compute_polytope_halfspaces
from numpy import exp,arange , arctan , sqrt, array
from scipy import signal
from math import cos,sin
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from subprocess import call
from matplotlib import cm
from shapely.geometry import Point
from shapely.geometry import LineString
from shapely.geometry.polygon import Polygon
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

# this is a parallelotope. 
#  ----5  -  7
#  |   /   / |
#  |  1 - 3  |
#  |  |   |  |
#  |  0 - 2  |            0 is the bottom left! meaning minimum x,y,x3 
#  |    \  \ |            7 is the up right! meaning maximum x,y,x3
#  |-----4 - 6 
# the links needed are:
# 0 -> 1 , 0 -> 2 , 0-> 4, 1-> 3 , 1->5, 2->3, 2->6, 3->7, 4->5, 4->6, 5->7, 6->7 
#


# some "global" variables
plot = True     # True for plotting
plot_volumes = True # True for step by step plotting of flow and estimation set parallelotopes
TRAN = True    # True for tran motion at step 1
plot_obstacles = False 


r = 6.6 /2.0    # 6.6(cm) radius of wheels in cm 
l = 13.2        # (cm) distance between the two weels 
if TRAN == True : 
    x1t = 12
    x2t = 13
else:
    x1t = random.randint(15,25)  # starting position for x1
    x2t = random.randint(15,25)  # starting position for x2
M = 0.5           # c boundary for M set
e = 0.2          # break condition for distance from target position
k1 = 10.8#10.8  # control gain for rotational motion
k2 = 1#0.05   # control gain for x1 in translational motion
k3 = 1#0.05   # control gain for x2 in translational motion

# some initialization of lists for printing
x1_set = []     # plot representative of x1 over time 
x2_set = []     # plot representative of x2 over time 
x3_set = []     # plot representative of x3 over time  
d_set = []      # plot d over time  
m_set = []      # plot m over time
df_set = []     # plot df over time (difference in angle)
angledf_set = []  # plot df over time (difference in angle) nikos method
temp_set =[]    # plot difference of two methods of df calculation over time 
cx = []         
cy = []
cz = []
hullx = []
hully = []
hullz = []
flow_volume = []
flow_parallelotope = []
estimation_volume = []
estimation_parallelotope = []
d1 = []
d2 = []
d3 = [] 
# some variables
run_time = 0.01      # discretization time 
discretization = 0.01
break_steps = 35    # break time! 
plotvar = 25         # define after how many steps the cube will be plotted! 
camera_steps = 10000   # define after how many steps we use the camera. Later on this maybe regarding the volume of the box
# static disturbances
d1_rot = [-0.0,0.0]#[-0.1,0.1]
d1_tran = [-0.1,0.1]#[-0.1,0.1]
d2_rot = [-0.0,0.0] #[-0.1,0.1]
d2_tran = [-0.1,0.1] #[-0.1,0.1]
d3_rot = [-0.1,0.1] #[-0.1,0.1]
d3_tran = [-1,1]
d1_camera = 0.00
d2_camera = 0.00
d3_camera = 0.00
# compatibility bounds
w1 = 0.0001  # static at the time being
w2 = 0.3   # 10% percentage of control action for angle

network_delay_distance = [304767.2202,220990.0076,169694.0005,135340.7176,110991.0751,93013.30647,79315.69879,68611.65572,60069.94168]
ap_locations = [[0,0],[25,25],[3,23],[23,5],[12,12]]

#TODO obstacles should be generated randomly. obstacles dont care about the orientation. should be columns in the 3D space. or plot only in 2d 
# obstacles are defined by two points according to: up left,  bottom right vertices of the obstacle (that is a parallelotope)
obstacles = [[[3,5],[4,4]],[[13,14],[14,12]]]
grid = [[[-1,0],[0,-1]],[[-1,25],[0,-1]],[[25,0],[26,-1]],[[25,26],[26,25]]] # make grid look as 4 obstacles
collision_threshold = 15 # threshold to assume that we are close to obstacles (in cms)
#TODO general TODOs
#TODO the hull lists are the same with estimation_parallelotope. Fix that. 
#TODO get_numbers should change if we want to use different run_time of flow and steps of flow.. now it is the same
#TODO when i use the camera plot trajectories and plot cube should change!!  
#TODO initial set should be an interval for x_3


def calculateslope(x1,x2):
    correct_angle = (np.arctan2(x2-x2t,x1-x1t))
    if  x1 >= x1t or (x1<x1t and x2<x2t): correct_angle += math.pi
    elif  x2 >= x2t : correct_angle = (np.arctan2(x2t-x2,x1t-x1))
    if correct_angle <0 : correct_angle = 2*math.pi - abs(correct_angle)
    return correct_angle

def rotational_noisy(est_set):
    x1 = [item[0] for item in est_set]
    x2 = [item[1] for item in est_set]
    x3 = [item[2] for item in est_set]
    xmin = min(x1)
    xmax = max(x1)
    ymin = min(x2)
    ymax = max(x2)
    x3min = min(x3)
    x3max = max(x3)

    u2 = cp.Variable()
    slopeopt = cp.Variable()
    x3cp = cp.Variable()
    
    slopes = []
    slopes.append(calculateslope(xmin,ymax))
    slopes.append(calculateslope(xmax,ymin))
    slopes.append(calculateslope(xmax,ymax))
    slopes.append(calculateslope(xmin,ymin))
    print ("slopes: ", slopes)
    slope_min = min(slopes)
    slope_max = max(slopes)
    if slope_min > slope_max: 
        temp = slope_max
        slope_max = slope_min
        slope_min = temp
    if slope_max > 5.5 and slope_min < 1:
        slope_max =  2*math.pi - slope_max  
    if slope_min > slope_max: 
        temp = slope_max
        slope_max = slope_min
        slope_min = temp

    prob = cp.Problem(cp.Minimize(cp.max(x3cp - slopeopt + u2 + d3_rot[0])  ),
    [               x3cp >= x3min , x3cp <= x3max, 
                    slopeopt >= slope_min, slopeopt <= slope_max, 
                    u2 >= -5 , u2 <= 10,
                    x3cp - slopeopt + u2 >= 0 ,
                    
    ])
    prob.solve()
    print("\nThe optimal value is", prob.value)
    print (x3cp.value, slopeopt.value, u2.value)
    x3 = x3cp.value + u2.value
    return u2.value, x3 
        
def translational_noisy(est_set):
    x1 = [item[0] for item in est_set]
    x2 = [item[1] for item in est_set]
    x3 = [item[2] for item in est_set]
    d1 = d1_tran
    d2 = d2_tran
    d3 = d3_tran
    xstar = [x1t,x2t]
    u = cp.Variable()

    prob = cp.Problem(cp.Minimize(cp.maximum(
    cp.norm2(cp.vstack([x1[0] + cos(x3[0])* u + d1[0] - xstar[0], x2[0] + sin(x3[0])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[0] + cos(x3[0])* u + d1[0] - xstar[0], x2[0] + sin(x3[0])*u + d2[1] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[0] + cos(x3[0])* u + d1[1] - xstar[0], x2[0] + sin(x3[0])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[0] + cos(x3[0])* u + d1[1] - xstar[0], x2[0] + sin(x3[0])*u + d2[1] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[0] + cos(x3[1])* u + d1[0] - xstar[0], x2[0] + sin(x3[1])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[0] + cos(x3[1])* u + d1[0] - xstar[0], x2[0] + sin(x3[1])*u + d2[1] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[0] + cos(x3[1])* u + d1[1] - xstar[0], x2[0] + sin(x3[1])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[0] + cos(x3[1])* u + d1[1] - xstar[0], x2[0] + sin(x3[1])*u + d2[1] -xstar[1]])) ,  
    cp.norm2(cp.vstack([x1[0] + cos(x3[0])* u + d1[0] - xstar[0], x2[1] + sin(x3[0])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[0] + cos(x3[0])* u + d1[0] - xstar[0], x2[1] + sin(x3[0])*u + d2[1] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[0] + cos(x3[0])* u + d1[1] - xstar[0], x2[1] + sin(x3[0])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[0] + cos(x3[0])* u + d1[1] - xstar[0], x2[1] + sin(x3[0])*u + d2[1] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[0] + cos(x3[1])* u + d1[0] - xstar[0], x2[1] + sin(x3[1])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[0] + cos(x3[1])* u + d1[0] - xstar[0], x2[1] + sin(x3[1])*u + d2[1] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[0] + cos(x3[1])* u + d1[1] - xstar[0], x2[1] + sin(x3[1])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[0] + cos(x3[1])* u + d1[1] - xstar[0], x2[1] + sin(x3[1])*u + d2[1] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[1] + cos(x3[0])* u + d1[0] - xstar[0], x2[0] + sin(x3[0])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[1] + cos(x3[0])* u + d1[0] - xstar[0], x2[0] + sin(x3[0])*u + d2[1] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[1] + cos(x3[0])* u + d1[1] - xstar[0], x2[0] + sin(x3[0])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[1] + cos(x3[0])* u + d1[1] - xstar[0], x2[0] + sin(x3[0])*u + d2[1] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[1] + cos(x3[1])* u + d1[0] - xstar[0], x2[0] + sin(x3[1])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[1] + cos(x3[1])* u + d1[0] - xstar[0], x2[0] + sin(x3[1])*u + d2[1] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[1] + cos(x3[1])* u + d1[1] - xstar[0], x2[0] + sin(x3[1])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[1] + cos(x3[1])* u + d1[1] - xstar[0], x2[0] + sin(x3[1])*u + d2[1] -xstar[1]])) ,  
    cp.norm2(cp.vstack([x1[1] + cos(x3[0])* u + d1[0] - xstar[0], x2[1] + sin(x3[0])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[1] + cos(x3[0])* u + d1[0] - xstar[0], x2[1] + sin(x3[0])*u + d2[1] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[1] + cos(x3[0])* u + d1[1] - xstar[0], x2[1] + sin(x3[0])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[1] + cos(x3[0])* u + d1[1] - xstar[0], x2[1] + sin(x3[0])*u + d2[1] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[1] + cos(x3[1])* u + d1[0] - xstar[0], x2[1] + sin(x3[1])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[1] + cos(x3[1])* u + d1[0] - xstar[0], x2[1] + sin(x3[1])*u + d2[1] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[1] + cos(x3[1])* u + d1[1] - xstar[0], x2[1] + sin(x3[1])*u + d2[0] -xstar[1]])) ,
    cp.norm2(cp.vstack([x1[1] + cos(x3[1])* u + d1[1] - xstar[0], x2[1] + sin(x3[1])*u + d2[1] -xstar[1]])) ,
                    )), 
                [u >= 0 , u<= 10])
    prob.solve()
    print ("Next u:", u.value)
    return u.value 

def rotational(x1,x2,x3,df):
    gain = r/l * k1 * df  * run_time
    x3 = x3 + gain 
    return x1,x2,x3,0,0,gain
    
def translational(x1,x2,x3):
    gainx1 = run_time *r * k2* math.sqrt((x1-x1t)**2+(x2-x2t)**2) 
    gainx2 = run_time *r * k3* math.sqrt((x1-x1t)**2+(x2-x2t)**2) 
    x1f = x1 +  gainx1 * math.cos(x3) 
    x2f = x2 +  gainx2 * math.sin(x3) 
    distance = (x1f-x1)**2 + (x2f-x2)**2
    return x1f,x2f,x3, gainx1, gainx2, 0, distance

def check_angle(x3):
    if x3 > 2*math.pi : 
        x3 = x3 - 2*math.pi
    if x3 <= 0:
        x3 = 2*math.pi + x3
    return x3 

def check_grid(x):
    if x > 250 :
        x = 250
    if x < 0:
        x = 0 
    return x 

def calculateangle(x1,x2,x3):
    xx = x1 + 1*math.cos(x3) 
    yy = x2 + 1*math.sin(x3)
	# x goes for j , y goes for i 	
	# A = (x,y) , B = (x1,y1) , C= (xtarget,ytarget)
	# find distance of the sides a_side = BC , b_side = AC , c_side = AB
    a_side = math.sqrt((xx-x1t)**2 + (yy-x2t)**2) 
    b_side = math.sqrt((x1-x1t)**2 + (x2-x2t)**2)
    c_side = math.sqrt((x1-xx)**2 + (x2-yy)**2)
    #print (a_side,b_side,c_side)
    #print ((a_side**2) + c_side**2 + b_side**2)
    #print (2*c_side*b_side)
    temp = (float(-(a_side**2) + c_side**2 + b_side**2)/(float(2*c_side*b_side))) 
	# Law of cosines a**2 = c**2 + b**2 - 2cbcosw , where w is the angle of A
    if temp > 1 : temp = 1
    w = math.acos(temp) 
    # We need also the slope to determine to whick way the Robot must rotate
    a = np.array([xx-x1,yy-x2])
    b = np.array([x1t-x1,x2t-x2])
    sign =  np.cross(a,b) # cross product of the two vectors
    #print (sign)
    #if sign != 0 :
    no =  w * sign/abs(sign)
    #else: 
	#    no = 0 
    return no

def plots(steps):
    t = np.linspace(0,steps,steps)
    #x1 and x2 together
    plt.plot(x1_set, x2_set, 'b.',label='translational motion over time')
    plt.xlabel('x1 (cm)', color='#1C2833')
    plt.ylabel('x2 (cm)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("./output/grid.png")
        #plt.show()
        plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax = plt.subplots()
    for est_set in estimation_parallelotope:
        listx = [item[0] for item in est_set]
        listy = [item[1] for item in est_set]
        height = max(listy) - min(listy)
        width = max(listx) -min(listx)
        ax.add_patch(Rectangle((min(listx), min(listx)), width, height,color ='green'))
    # ax.xlabel('x1 (cm)', color='#1C2833')
    # ax.ylabel('x2 (cm)', color='#1C2833')
    # plt.legend(loc='upper left')
    #display plot
    plt.xlim([5, 15])
    plt.ylim([5, 15])
    plt.show()

    # # x1
    # plt.plot(t, x1_set, 'g.',label='x1 over time')
    # plt.xlabel('steps', color='#1C2833')
    # plt.ylabel('x1 (cm)', color='#1C2833')
    # plt.legend(loc='upper left')
    # plt.grid()
    # if plot:
    #     plt.savefig("./output/x1.png")
    #     #plt.show()
    #     plt.close()

    # # x2
    # plt.plot(t, x2_set, 'b.',label='x2 over time')
    # plt.xlabel('steps', color='#1C2833')
    # plt.ylabel('x2 (cm)', color='#1C2833')
    # plt.legend(loc='upper left')
    # plt.grid()
    # if plot:
    #     plt.savefig("./output/x2.png")
    #     #plt.show()
    #     plt.close()

    # # x3
    # plt.plot(t, x3_set, 'r.',label='x3 over time')
    # plt.xlabel('steps', color='#1C2833')
    # plt.ylabel('x3 (degrees)', color='#1C2833')
    # plt.legend(loc='upper left')
    # plt.grid()
    # if plot:
    #     plt.savefig("./output/x3.png")
    #     #plt.show()
    #     plt.close()

    # m
    # plt.plot(t, m_set, 'b.',label='m over time/ limit=%d' %M)
    # plt.axhline(M,color='red')
    # plt.xlabel('steps', color='#1C2833')
    # plt.ylabel('m ', color='#1C2833')
    # plt.legend(loc='upper left')
    # plt.grid()
    # if plot:
    #     plt.savefig("./output/m.png")
    #     #plt.show()
    #     plt.close()

    # d

    plt.plot(t, d_set, 'b.',label='d over time/ limit=%d' %e)
    plt.axhline(e,color='red')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('d (cm)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("./output/d.png")
        #plt.show()
        plt.close()

    # # angle dif 
    # plt.plot(t, angledf_set, 'r.',label='angledif over time')
    # plt.xlabel('steps', color='#1C2833')
    # plt.ylabel('angledif  (degrees)', color='#1C2833')
    # plt.legend(loc='upper left')
    # plt.grid()
    # if plot:
    #     plt.savefig("./output/angledif.png")
    #     #plt.show()
    #     plt.close()

    # # df
    # plt.plot(t, df_set, 'r.',label='angle needed to rotate in order to look target')
    # plt.xlabel('steps', color='#1C2833')
    # plt.ylabel('angle needed to rotate in order to look target (degrees)', color='#1C2833')
    # plt.legend(loc='upper left')
    # plt.grid()
    # if plot:
    #     plt.savefig("./output/df.png")
    #     #plt.show()
    #     plt.close()

    # # difference between two methods
    # plt.plot(t, temp_set, 'r.',label='difference between two methods')
    # plt.xlabel('steps', color='#1C2833')
    # plt.ylabel('difference between two methods (degrees)', color='#1C2833')
    # plt.legend(loc='upper left')
    # plt.grid()
    # if plot:
    #     plt.savefig("./output/methods.png")
    #     #plt.show()
    #     plt.close()

    t = np.linspace(0,steps+1,steps+1)
    # flow volume, estimation_volume
    plt.plot(t, flow_volume, 'r.',label='volumes of flow set parallelotope')
    plt.plot(t, estimation_volume, 'b.',label='volumes of estimation set parallelotope cm^3')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('volumes of flow , estimation set at each step', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("./output/volumes.png")
        #plt.show()
        plt.close()

def plot_cube():
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_zlabel('x3')

    for i in range (0,8):
        printx = []
        printy = []
        printx3 = []
        for x in range (0 , len(hullx),plotvar): 
            printx.append(hullx[x][i])      # ta x ana 10 kai ana shmeio tou kubou 
            printy.append(hully[x][i])
            printx3.append(hullz[x][i])
        ax.scatter(printx, printy, printx3, c='b',s=20)
        ax.plot(printx, printy, printx3, color= 'r')
    
    for o in range (0 , len(hullx),plotvar): 
        for i in range (0,7,2): 
            printx = []
            printy = []
            printx3 = []
            printx.append(hullx[o][i])      
            printx.append(hullx[o][i+1])
            printy.append(hully[o][i])
            printy.append(hully[o][i+1])
            printx3.append(hullz[o][i])
            printx3.append(hullz[o][i+1])
            #sthn plot pane ta x values duo shmeiwn kai meta ta y values 2 shmeiwn
            ax.plot(printx, printy, printx3, color=   'g')

        for i in range (0,6):
            if (i!=2 and i!=3): 
                printx = []
                printy = []
                printx3 = []
                printx.append(hullx[o][i])      
                printx.append(hullx[o][i+2])
                printy.append(hully[o][i])
                printy.append(hully[o][i+2])
                printx3.append(hullz[o][i])
                printx3.append(hullz[o][i+2])
                #sthn plot pane ta x values duo shmeiwn kai meta ta y values 2 shmeiwn
                ax.plot(printx, printy, printx3, color=   'g')
        
        for i in range (0,4):
            printx = []
            printy = []
            printx3 = []
            printx.append(hullx[o][i])      
            printx.append(hullx[o][i+4])
            printy.append(hully[o][i])
            printy.append(hully[o][i+4])
            printx3.append(hullz[o][i])
            printx3.append(hullz[o][i+4])
            #sthn plot pane ta x values duo shmeiwn kai meta ta y values 2 shmeiwn
            ax.plot(printx, printy, printx3, color=   'g')
    
    # center of parallelotope 
    cxlist = []
    cylist = []
    cx3list = []
    for x in range (0 , len(hullx),plotvar): 
        cxlist.append(cx[x])      # ta x ana 10 kai ana shmeio tou kubou 
        cylist.append(cy[x])
        cx3list.append(cz[x])
    ax.scatter3D(cxlist, cylist, cx3list)
    plt.show()  

def get_volumes(xminstar,xmaxstar,yminstar,ymaxstar,x3minstar,x3maxstar,x_under,x_over,y_under,y_over,x3_under,x3_over,u1, u2, u3, d1, d2, d3):
    estimation_volume.append((xmaxstar-xminstar)*(ymaxstar-yminstar)*(x3maxstar-x3minstar))
    flow_volume.append((x_over-x_under)*(y_over-y_under)*(x3_over-x3_under))
    estimation_parallelotope.append(list(itertools.product([xminstar,xmaxstar],[yminstar,ymaxstar],[x3minstar,x3maxstar])))
    flow_parallelotope.append(list(itertools.product([x_under,x_over],[y_under,y_over],[x3_under,x3_over])))

    if (plot_volumes == True):
        fig_vol = plt.figure()
        vol = fig_vol.gca(projection='3d')
        vol.set_xlabel('x1')
        vol.set_ylabel('x2')
        vol.set_zlabel('x3')
        for i in range (0,7,2): 
            vol1x = []
            vol1y = []
            vol1z = []
            vol2x = []
            vol2y = []
            vol2z = []
            vol3x = []
            vol3y = []
            vol3z = []
            vol1x.append(estimation_parallelotope[-1][i][0])      
            vol1x.append(estimation_parallelotope[-1][i+1][0])
            vol1y.append(estimation_parallelotope[-1][i][1])
            vol1y.append(estimation_parallelotope[-1][i+1][1])
            vol1z.append(estimation_parallelotope[-1][i][2])
            vol1z.append(estimation_parallelotope[-1][i+1][2])
            if i == 0 :
                vol.plot(vol1x, vol1y, vol1z, color= 'g', label = "New Estimation Set")
            else: 
                vol.plot(vol1x, vol1y, vol1z, color= 'g')
            vol2x.append(flow_parallelotope[-1][i][0])      
            vol2x.append(flow_parallelotope[-1][i+1][0])
            vol2y.append(flow_parallelotope[-1][i][1])
            vol2y.append(flow_parallelotope[-1][i+1][1])
            vol2z.append(flow_parallelotope[-1][i][2])
            vol2z.append(flow_parallelotope[-1][i+1][2])
            if i == 0 :
                vol.plot(vol2x, vol2y, vol2z, linestyle=':' , color='r',  label = "TM over-approximation")
            else: 
                vol.plot(vol2x, vol2y, vol2z, linestyle=':' , color='r')
            vol3x.append(estimation_parallelotope[-2][i][0])      
            vol3x.append(estimation_parallelotope[-2][i+1][0])
            vol3y.append(estimation_parallelotope[-2][i][1])
            vol3y.append(estimation_parallelotope[-2][i+1][1])
            vol3z.append(estimation_parallelotope[-2][i][2])
            vol3z.append(estimation_parallelotope[-2][i+1][2])
            if i == 0 :
                vol.plot(vol3x, vol3y, vol3z, linestyle=':' ,color= 'b', label = "Current Estimation Set")
            else: 
                vol.plot(vol3x, vol3y, vol3z, linestyle=':' ,color= 'b')

        for i in range (0,6):
            if (i!=2 and i!=3): 
                vol1x = []
                vol1y = []
                vol1z = []
                vol2x = []
                vol2y = []
                vol2z = []
                vol3x = []
                vol3y = []
                vol3z = []
                vol1x.append(estimation_parallelotope[-1][i][0])      
                vol1x.append(estimation_parallelotope[-1][i+2][0])
                vol1y.append(estimation_parallelotope[-1][i][1])
                vol1y.append(estimation_parallelotope[-1][i+2][1])
                vol1z.append(estimation_parallelotope[-1][i][2])
                vol1z.append(estimation_parallelotope[-1][i+2][2])
                vol.plot(vol1x, vol1y, vol1z, color=   'g')
                vol2x.append(flow_parallelotope[-1][i][0])      
                vol2x.append(flow_parallelotope[-1][i+2][0])
                vol2y.append(flow_parallelotope[-1][i][1])
                vol2y.append(flow_parallelotope[-1][i+2][1])
                vol2z.append(flow_parallelotope[-1][i][2])
                vol2z.append(flow_parallelotope[-1][i+2][2])
                vol.plot(vol2x, vol2y, vol2z, linestyle=':' ,color=   'r')
                vol3x.append(estimation_parallelotope[-2][i][0])      
                vol3x.append(estimation_parallelotope[-2][i+2][0])
                vol3y.append(estimation_parallelotope[-2][i][1])
                vol3y.append(estimation_parallelotope[-2][i+2][1])
                vol3z.append(estimation_parallelotope[-2][i][2])
                vol3z.append(estimation_parallelotope[-2][i+2][2])
                vol.plot(vol3x, vol3y, vol3z, linestyle=':' ,color= 'b')
        for i in range (0,4):
            vol1x = []
            vol1y = []
            vol1z = []
            vol2x = []
            vol2y = []
            vol2z = []
            vol3x = []
            vol3y = []
            vol3z = []
            vol1x.append(estimation_parallelotope[-1][i][0])      
            vol1x.append(estimation_parallelotope[-1][i+4][0])
            vol1y.append(estimation_parallelotope[-1][i][1])
            vol1y.append(estimation_parallelotope[-1][i+4][1])
            vol1z.append(estimation_parallelotope[-1][i][2])
            vol1z.append(estimation_parallelotope[-1][i+4][2])
            vol.plot(vol1x, vol1y, vol1z, color=   'g')
            vol2x.append(flow_parallelotope[-1][i][0])      
            vol2x.append(flow_parallelotope[-1][i+4][0])
            vol2y.append(flow_parallelotope[-1][i][1])
            vol2y.append(flow_parallelotope[-1][i+4][1])
            vol2z.append(flow_parallelotope[-1][i][2])
            vol2z.append(flow_parallelotope[-1][i+4][2])
            vol.plot(vol2x, vol2y, vol2z, linestyle=':' ,color= 'r')
            vol3x.append(estimation_parallelotope[-2][i][0])      
            vol3x.append(estimation_parallelotope[-2][i+4][0])
            vol3y.append(estimation_parallelotope[-2][i][1])
            vol3y.append(estimation_parallelotope[-2][i+4][1])
            vol3z.append(estimation_parallelotope[-2][i][2])
            vol3z.append(estimation_parallelotope[-2][i+4][2])
            vol.plot(vol3x, vol3y, vol3z, linestyle=':' ,color= 'b')
        # plot trajectories from previous estimation set vertices given the control action to validate that exist inside the new estimation set
        # for x in estimation_parallelotope[-2]:
        #     for i in range (2):
        #         for j in range (2):
        #             for k in range (2):
        #                 x_trajectory = []
        #                 y_trajectory = []
        #                 x3_trajectory = []
        #                 #TODO change this and use the function of control actions.. Crucial!
        #                 if (u3==0):
        #                     x_trajectory = [x[0],(x[0]+u1*math.cos(x[2]))+run_time *d1[i]]
        #                     y_trajectory = [x[1],(x[1]+u2*math.sin(x[2]))+run_time *d2[j]]
        #                     x3_trajectory = [x[2],(x[2]+run_time *d3[k])]
        #                 else :
        #                     x_trajectory = [x[0],(x[0]+run_time *d1[i])]
        #                     y_trajectory = [x[1],(x[1]+run_time *d2[j])]
        #                     #TODO x3_trajectory should have +d3
        #                     x3_trajectory = [x[2],(x[2]+u3)]
        #                 point = Point(x_trajectory[1],y_trajectory[1],x3_trajectory[1])
        #                 polygon = Polygon(flow_parallelotope[-1])
        #                 if (polygon.exterior.distance(point)) > 1e-2 : 
        #                     print ("Danger possible trajectory outside of flow!")
        #                     print (polygon.exterior.distance(point))
        #                 vol.plot(x_trajectory,y_trajectory,x3_trajectory,linestyle='--' ,color= 'y')
        
        if plot_obstacles:
            for item in obstacles:
                maxx = item[1][0]
                minx = item[0][0]
                maxy = item[1][1]
                miny = item[0][1]
                minx3 = x3minstar
                maxx3 = x3maxstar
                obstacle_parallelotope = []
                obstacle_parallelotope.append([minx,miny,minx3])
                obstacle_parallelotope.append([minx,miny,maxx3]) 
                obstacle_parallelotope.append([minx,maxy,minx3]) 
                obstacle_parallelotope.append([minx,maxy,maxx3]) 
                obstacle_parallelotope.append([maxx,miny,minx3]) 
                obstacle_parallelotope.append([maxx,miny,maxx3]) 
                obstacle_parallelotope.append([maxx,maxy,minx3]) 
                obstacle_parallelotope.append([maxx,maxy,maxx3]) 
                for i in range (0,7,2): 
                    volobx = []
                    voloby = []
                    volobz = []
                    volobx.append(obstacle_parallelotope[i][0])      
                    volobx.append(obstacle_parallelotope[i+1][0])
                    voloby.append(obstacle_parallelotope[i][1])
                    voloby.append(obstacle_parallelotope[i+1][1])
                    volobz.append(obstacle_parallelotope[i][2])
                    volobz.append(obstacle_parallelotope[i+1][2])
                    if i == 0 and obstacles.index(item)==0:
                        vol.plot(volobx, voloby, volobz, linestyle='--',color= 'k', label = "Obstacles")
                    else: 
                        vol.plot(volobx, voloby, volobz, linestyle='--',color= 'k')
                for i in range (0,6):
                    if (i!=2 and i!=3): 
                        volobx = []
                        voloby = []
                        volobz = []
                        volobx.append(obstacle_parallelotope[i][0])      
                        volobx.append(obstacle_parallelotope[i+2][0])
                        voloby.append(obstacle_parallelotope[i][1])
                        voloby.append(obstacle_parallelotope[i+2][1])
                        volobz.append(obstacle_parallelotope[i][2])
                        volobz.append(obstacle_parallelotope[i+2][2])
                        vol.plot(volobx, voloby, volobz, linestyle='--',color= 'k')
                for i in range (0,4):
                    volobx = []
                    voloby = []
                    volobz = []
                    volobx.append(obstacle_parallelotope[i][0])      
                    volobx.append(obstacle_parallelotope[i+4][0])
                    voloby.append(obstacle_parallelotope[i][1])
                    voloby.append(obstacle_parallelotope[i+4][1])
                    volobz.append(obstacle_parallelotope[i][2])
                    volobz.append(obstacle_parallelotope[i+4][2])
                    vol.plot(volobx, voloby, volobz, linestyle='--',color= 'k')

        # save figure
        plt.legend(loc="upper left")
        i = str(len(flow_parallelotope)-1)
        plt.savefig("./output/cube"+i+".png")     
        #plt.show()  
        plt.close()

def lists_renew(x1,x2,x3,distance,m,df,angle):
    x1_set.append(x1)
    x2_set.append(x2)
    x3_set.append(math.degrees(x3))
    d_set.append(distance)
    m_set.append(m)
    df_set.append(math.degrees(df))
    angledf_set.append(math.degrees(angle))
    temp_set.append(math.degrees(df)-math.degrees(angle))

def build_model(time, d1, d2, d3, u1, u2, u3, jumptime, initial_set, mode,output,name):
    with open('replan.model', 'w') as the_file:
        the_file.write('hybrid reachability\n''{\n''  state var x,y,x3,t \n'' setting\n''{\n'
        'fixed steps '+str(run_time)+'\n''time '+str(discretization)+'\n' # time  
        'remainder estimation 1e-2\n''identity precondition\n''gnuplot grid 5 '+ output+'\n'
        'fixed orders 8\n''cutoff 1e-15\n''precision 2000\n''output '+name+'\n''max jumps 1\n''print on\n''}\n')
        the_file.write('modes\n'
        '{\n''tran\n''{ \n''nonpoly ode\n''{\n'"t' = 1\n")
        the_file.write("x3' ="+ str(d3) + "\n" # to d 1
        "x' = "+str(u1)+"*cos(x3) +"+str(d1)+ "\n"   # to u 1
        "y' = "+str(u2)+"*sin(x3) +"+str(d2)+ "\n""}\n") # to u 2
        the_file.write("inv\n""{\n""t <= 1""\n""t<=1""\n") 
        # for item in invtran:
        #     the_file.write(item+"\n")
        the_file.write("\n"   
        "}\n""}\n"
        "rot\n""{\n""nonpoly ode\n""{\n"
        "x3' = "+str(u3)+"+"+str(d3)+"\n"  # to u3
        "x' = "+str(d1)+"\n""y' = "+str(d2)+"\n" # to d1
        "t' = 1\n""}\n""inv \n""{\n" # to d2
        "" "t <= 1""\n""t<=1""\n")
        # for item in invrot:
        #     the_file.write(item+"\n")
        the_file.write(" \n"
        "}\n""}\n""}\n"
        "jumps\n""{\n""tran -> rot\n"
        "guard { t = "+str(jumptime)+"  }\n" # edw einai h fasoula gia to jump
        "reset { }\n""parallelotope aggregation { }\n""rot -> tran \n""guard {x = 5}\n""reset{}#x' := x - 4.9 }\n""parallelotope aggregation { }\n""}\n"
        "init\n""{\n"""+mode +"\n""{\n"
        "x in "+initial_set[0]+"\n" # in "+initial_set[0]+"\n"  # add initial set constraints 
        "y in "+initial_set[1]+"\n"  # add initial set constraints 
        "x3 in "+initial_set[2]+"\n"   # add initial set constraints  
        "}\n""}\n""}\n"
    )
    ["t <= 1"] , ["t<=1"]

def get_numbers():
    octagon_vertices=[]
    with open("octagon.txt") as f:
        for line in f.readlines()[10:-3]:    # TODO this is very manual. we should target the last lines of the file if we choose longer run_time
            if line == '\n': continue; 
            x=[]
            index = 0 
            b=0
            for item in line.split():
                if index == 0 :
                    a = item
                else: 
                    b = float(item)
                index +=1
            if b == float(run_time):
                x.append(float(a))
                octagon_vertices.append(x)
    return min(octagon_vertices)[0], max(octagon_vertices)[0]

def get_initialset(est_set):
    listx = [item[0] for item in est_set]
    listy = [item[1] for item in est_set]
    listx3 = [item[2] for item in est_set]
    initial_set = ["["+str(min(listx))+","+str(max(listx))+"]","["+str(min(listy))+","+str(max(listy))+"]","["+str(min(listx3))+","+str(max(listx3))+"]",]
    # print ("Initial set: ")
    # print (initial_set)
    # input("Press Enter to continue...")
    return initial_set, min(listx),max(listx), min(listy),max(listy), min(listx3),max(listx3)

def check_feasibility(est_set):
    listx3 = [item[2] for item in est_set]
    minx3 = min(listx3)
    maxx3 = max(listx3)
    no_borders = obstacles + grid
    collision = False
    distance_obstacle = 15000
    vertice_danger = None
    closest_point = None 
    for obstacle in no_borders:
        maxx = obstacle[1][0]
        minx = obstacle[0][0]
        maxy = obstacle[1][1]
        miny = obstacle[0][1]
        obstacle_parallelotope = []
        obstacle_parallelotope.append([minx,miny,minx3])
        obstacle_parallelotope.append([minx,miny,maxx3]) 
        obstacle_parallelotope.append([minx,maxy,minx3]) 
        obstacle_parallelotope.append([minx,maxy,maxx3]) 
        obstacle_parallelotope.append([maxx,miny,minx3]) 
        obstacle_parallelotope.append([maxx,miny,maxx3]) 
        obstacle_parallelotope.append([maxx,maxy,minx3]) 
        obstacle_parallelotope.append([maxx,maxy,maxx3])
        for vertice in est_set:
            point = Point(vertice[0],vertice[1],vertice[2])
            polygon = Polygon(obstacle_parallelotope)
            if (polygon.exterior.distance(point)) < collision_threshold : 
                if distance_obstacle > polygon.exterior.distance(point): 
                    distance_obstacle = polygon.exterior.distance(point)
                    #vertice_danger = est_set.index(vertice)
                    collision = True
                    closest_obstacle = no_borders.index(obstacle)  # TODO useless
                    # Define half space representation and solve a simple optimization problem with infinity norm
                    #TODO only do that once! 
                    vertices = map(array, obstacle_parallelotope)
                    A, b = compute_polytope_halfspaces(vertices)
                    vertices = map(array, est_set)
                    C, d = compute_polytope_halfspaces(vertices)
                    x = cp.Variable(len(A[0]))
                    y = cp.Variable(len(C[0]))
                    prob1 = cp.Problem(cp.Minimize(cp.norm_inf(x-y)),
                                        [A*x <= b ,
                                        C*y <= d])
                    try:
                        distance_obstacle = prob1.solve()
                    except Exception as e:
                        print(e)
                    closest_point = y.value 

    if collision:
        print ("Danger close to ostacle!")
        print ("Distance to obstacle: " "{:10.2f}".format(distance_obstacle))
        print ("Point (x,y,x3) of estimation set: ",closest_point)
    return collision, closest_point

def get_tm_intervals(mode): 
    with open("flopipes.txt") as f:
        for line in f.readlines()[10:-3]:    # TODO this is very manual. we should target the last lines of the file if we choose longer run_time
            if 'x = ' in line:
                x_calc = line 
            if 'y = ' in line:
                y_calc = line
            if 'x3 = ' in line:
                x3_calc = line 
    if (mode == "tran"):
        x_flow_min , x_flow_max = get_nonlinear_flowpipe("x",x_calc)
        y_flow_min , y_flow_max = get_nonlinear_flowpipe("y",y_calc)
        x3_flow_min , x3_flow_max = get_linear_flowpipe("x3",x3_calc)

    else: 
        x_flow_min , x_flow_max = get_linear_flowpipe("x",x_calc)
        y_flow_min , y_flow_max = get_linear_flowpipe("y",y_calc)
        x3_flow_min , x3_flow_max = get_linear_flowpipe("x3",x3_calc)
    return x_flow_min,x_flow_max,y_flow_min,y_flow_max,x3_flow_min,x3_flow_max

def get_nonlinear_flowpipe (state_var,state_calc): 
    local_t = run_time
    local_var_1 = np.linspace(-1 , 1, 100)
    local_var_2 = np.linspace(-1 , 1, 100)
    local_var_3 = np.linspace(-1 , 1, 100)
    local_var_4 = np.linspace(-1 , 1, 100)
    #min 
    state_safe = state_calc
    state_calc = state_calc.replace("[","min(")   
    state_calc = state_calc.replace("]",")")    
    state_calc = state_calc.replace("^","**")   
    state_calc = state_calc.replace(""+state_var+" = ","")   
    state_calc = state_calc.replace("\n","") 
    state_calc_min = ""
    for item in state_calc.split(' ') :
        if (item.find("min")) != -1 :
            if item[4]=='-': 
                item = item.replace("min","max")
        state_calc_min += (item)
    location = (state_calc_min.find("local_var_"))
    #print (location)
    if location == -1:   # catch the case where it is not actually a non_linear flowpiper after all! dunnoy why
        return get_linear_flowpipe(state_var,state_calc)
    try:
        state_a_min = min(eval(state_calc_min[:location+11]))
    except TypeError:
        state_a_min = eval(state_calc_min[:location+11])
    try:
        state_b_min = min(eval(state_calc_min[location+12:]))
    except TypeError:
        state_b_min = eval(state_calc_min[location+12:])

    state_out_min =  (state_a_min +state_b_min)
    #print (state_out_min) 
    #max
    state_calc = state_safe
    state_calc = state_calc.replace("[","max(")   
    state_calc = state_calc.replace("]",")")    
    state_calc = state_calc.replace("^","**")   
    state_calc = state_calc.replace(""+state_var+" = ","") 
    state_calc = state_calc.replace("\n","") 
    state_calc_min = ""
    for item in state_calc.split(' ') :
        if (item.find("max")) != -1 : 
            if item[4]=='-': 
                item = item.replace("max","min")
        state_calc_min += (item)
    location = (state_calc_min.find("local_var_"))
    try:
         state_a_max = max(eval(state_calc_min[:location+11]))
    except TypeError:
         state_a_max = eval(state_calc_min[:location+11])
    try:
        state_b_max = max(eval(state_calc_min[location+12:]))
    except TypeError:
        state_b_max = eval(state_calc_min[location+12:])
    state_out_max =  (state_a_max +state_b_max)
    #print (state_out_max) 
    return state_out_min, state_out_max

def get_linear_flowpipe(state_var, state_calc):
    # print (state_var)
    # print (state_calc)
    local_t = run_time
    local_var_1 = np.linspace(-1 , 1, 100)
    local_var_2 = np.linspace(-1 , 1, 100)
    local_var_3 = np.linspace(-1 , 1, 100)
    local_var_4 = np.linspace(-1 , 1, 100)
    #min
    state_safe = state_calc
    state_calc = state_calc.replace("[","max(")   
    state_calc = state_calc.replace("]",")")    
    state_calc = state_calc.replace("^","**")   
    state_calc = state_calc.replace(state_var+" = ","")  
    state_calc = state_calc.replace("\n","")
    try:
        state_out_max = max(eval(state_calc))
    except TypeError:
        state_out_max = eval(state_calc)
    #print (state_calc)
    #max
    state_calc = state_safe
    state_calc = state_calc.replace("[","min(")   
    state_calc = state_calc.replace("]",")")    
    state_calc = state_calc.replace("^","**")   
    state_calc = state_calc.replace(state_var+" = ","")  
    state_calc = state_calc.replace("\n","") 
    try:
        state_out_min = min(eval(state_calc))
    except TypeError:
        state_out_min = eval(state_calc)
    #print (x_out_min)
    #print (x_out_max)
    return state_out_min, state_out_max

def compatibility(time, d1, d2, d3, u1, u2, u3, jumptime, initial_set, mode,xmin,xmax,ymin,ymax,x3min,x3max,distance,x3,x_under,x_over,y_under,y_over,x3_under,x3_over):
    measur1 = distance
    y = distance 

    wtran =  w1 * distance
    wangle = w2 * 5
    #wangle = random.uniform(0.1,0.3)
    #let's find x1,x2 intervals. the flow results cut with compatibility equations chosen. 
    yrootmin, yrootmax = optimize_quadratic(y_under,y_over, ymin,ymax,measur1,wtran)   # I want the minimum and maximum value of under root! check overleaf
    xrootmin, xrootmax = optimize_quadratic(x_under,x_over, xmin,xmax,measur1,wtran)   
    
    if xmin >= x_under :    #TODO edw thelei allagi logika. kati den paei kala! 
        #print ("negative for x's")
        xminstar = xmin - math.sqrt(y-wtran -yrootmin) 
        if xminstar < x_under : xminstar = x_under
        xmaxstar = xmax - math.sqrt(y+wtran -yrootmax) 
        if xmaxstar > x_over : 
            xmaxstar = x_over 
        #if xmaxstar < xminstar : xmaxstar = xminstar

    else : 
        if y-wtran -yrootmax < 0 : xminstar = xmin
        else:
            xminstar = xmin + math.sqrt(y-wtran -yrootmax) #max(x_under, (xmin))   # wrong if negative motion!! 
        if xminstar < x_under : xminstar = x_under
        #if x_under < xmax: xrootmin = 0 
        xmaxstar = xmax + math.sqrt(y+wtran -yrootmin) #min(x_over, (xmax +math.sqrt(measur1+wtran))
        if xmaxstar > x_over : 
            #print ("panw orio", xmaxstar, yrootmin, math.sqrt(y+wtran -yrootmin) )
            xmaxstar = x_over 

    if ymin >= y_under :
        #print ("negative for y's")
        if y-wtran -xrootmin < 0 : yminstar = ymin
        yminstar = ymin - math.sqrt(y-wtran -xrootmin) 
        if yminstar < y_under : yminstar = y_under
        ymaxstar = ymax - math.sqrt(y+wtran -xrootmax)  
        if ymaxstar > y_over : ymaxstar = y_over 
        #if ymaxstar < yminstar : ymaxstar = yminstar

    else: 
        print (y , wtran, xrootmax)
        yminstar = ymin + math.sqrt(y-wtran -xrootmax) 
        if yminstar < y_under : yminstar = y_under
        #if y_under < ymax: yrootmin = 0 
        ymaxstar = ymax + math.sqrt(y+wtran -xrootmin)  #min(y_over, (ymax +math.sqrt(measur1+wtran)))
        if ymaxstar > y_over : ymaxstar = y_over 
       
    # find x3 interval! 
    x3minstar, x3maxstar = optimization_problem(x3_under,x3_over,x3min,x3max,u3,wangle)    
    print ("Optimization", xminstar,xmaxstar,yminstar,ymaxstar,x3minstar,x3maxstar)

    xvalue,yvalue,x3value,cube = calculateEstimationSet(xminstar,xmaxstar,yminstar,ymaxstar,x3minstar,x3maxstar)
    get_volumes(xminstar,xmaxstar,yminstar,ymaxstar,x3minstar,x3maxstar,x_under,x_over,y_under,y_over,x3_under,x3_over,u1, u2, u3, d1, d2, d3)
    return xvalue,yvalue,x3value,cube

def run_flow(time, d1, d2, d3, u1, u2, u3, jumptime, initial_set, mode):
    output = "x,t"
    build_model(time, d1, d2, d3, u1/run_time, u2/run_time, u3/run_time, jumptime , initial_set, mode , output, "x") 
    #build_model(time, d1, d2, d3, u1, u2, u3, jumptime , initial_set, mode , output, "x") 

    rc = call(["./run.sh","x"])
    # #verticesx = get_numbers()
    # x_under, x_over = get_numbers()
    # #input("Press Enter to continue...")
    # output = "y,t"
    # build_model(time, d1, d2, d3, u1/run_time, u2/run_time, u3/run_time,  jumptime, initial_set, mode , output, "y") 
    # rc = call(["./run.sh","y"])
    # #verticesy = get_numbers()
    # y_under, y_over = get_numbers()
    # #input("Press Enter to continue...")
    # output = "x3,t"
    # build_model(time, d1, d2, d3, u1/run_time, u2/run_time, u3/run_time,jumptime, initial_set, mode , output, "x3") 
    # rc = call(["./run.sh","x3"])
    # #verticesx3 = get_numbers()
    # x3_under, x3_over = get_numbers()
    ###########
    x_under,x_over,y_under,y_over,x3_under,x3_over = get_tm_intervals(mode)
    print (mode)
    print ("Gains: ",u1,u2,u3)
    #input("Press Enter to continue...")
    print ("Flow:",x_under,x_over,y_under,y_over,x3_under,x3_over)
    return x_under,x_over,y_under,y_over,x3_under,x3_over

def optimization_problem(v3min,v3max,z3min,z3max,y,w):
    v3 = cp.Variable()
    z3 = cp.Variable()
    prob1 = cp.Problem(cp.Minimize(v3),
                        [v3 - z3   >= y-w , 
                        v3 - z3   <= y+w , 
                        z3 >= z3min , z3 <= z3max, 
                        v3 >= v3min , v3 <= v3max ])
    
    try:
        minimum = prob1.solve()
    except Exception as e:
        print(e)
    v3 = cp.Variable()
    z3 = cp.Variable()
    prob1 = cp.Problem(cp.Maximize(v3),
                        [v3 - z3   >= y-w , 
                        v3 - z3   <= y+w , 
                        z3 >= z3min , z3 <= z3max, 
                        v3 >= v3min , v3 <= v3max ])
    
    try:
        maximum = prob1.solve()
    except Exception as e:
        print(e)
    return minimum , maximum

def optimize_quadratic(vmin,vmax,zmin,zmax,y,w):
    #print (vmin,vmax,zmin,zmax,y,w)
    if zmax >= vmax : w = -w 
    v = cp.Variable()
    z = cp.Variable()
    prob1 = cp.Problem(cp.Minimize(cp.square(v-z)),
                        [cp.square(v-z)   <= y+w, 
                        z >= zmin , z <= zmax, 
                        v >= vmin , v <= vmax ])
    
    try:
        minimum = prob1.solve()
    except Exception as e:
        print(e)

    v = cp.Variable()
    z = cp.Variable()
    prob1 = cp.Problem(cp.Maximize(v-z),
                        [cp.square(v-z)   <= y-w, 
                        z >= zmin , z <= zmax, 
                        v >= vmin , v <= vmax ])
    if zmax >= vmax : 
        #print ("going negative")
        #going to negative of this state variable!! When we want to maximize a square and we have negative values. lets just minimize the inside! 
        prob1 = cp.Problem(cp.Minimize(v-z),   # minimize a convex is not possible so we do a trick 
                        [cp.square(v-z)   <= y-w, 
                        z >= zmin , z <= zmax, 
                        v >= vmin , v <= vmax ])
    try:
        maximum = prob1.solve()
    except Exception as e:
        print(e)   
    maximum = maximum**2  
    #print (minimum , maximum)
    if maximum >= y - w : 
        maximum = y - w    # needed for round of dcp problem!!  
    if minimum >= y + w : 
        minimum = 0        # TODO needs reasoning!! 
    #print (minimum , maximum)
    return minimum , maximum 

def calculateEstimationSet(xmin,xmax,ymin,ymax,x3min,x3max):  
    vertices = list(itertools.product([xmin,xmax],[ymin,ymax],[x3min,x3max]))
    listx = [item[0] for item in vertices]
    listy = [item[1] for item in vertices]
    listx3 = [item[2] for item in vertices]
    cx.append((xmin+xmax)/2)
    cy.append((ymin+ymax)/2) 
    cz.append((x3min+x3max)/2)
    hullx.append(listx)
    hully.append(listy)
    hullz.append(listx3)  
    return cx[-1],cy[-1],cz[-1], vertices

def calculateConvexHull(est_set,vertices):  
    #Instead of Convex Hull we get the center for each state variable.. 
    temp = []
    #start_time = time.time()
    for item in vertices: 
        temp.append([check_grid(item[0]),check_grid(item[1]),item[2]])
    #print ("vertices")
    #print (vertices)
    try:
        points = np.vstack(temp)
    except ValueError:
        print ("Value Error")
        print (temp)
        input("Press Enter to continue...")
    #print ("Time for vstack", time.time()-start_time)
    #start_time = time.time()
    hull = ConvexHull(points)
    #print ("Number of points :", len(points))
    #print ("Time for convex hull", time.time()-start_time)
    #Get centoid
    cx.append(np.mean(hull.points[hull.vertices,0]))
    cy.append(np.mean(hull.points[hull.vertices,1]))
    cz.append(np.mean(hull.points[hull.vertices,2]))
    tempx = np.mean(hull.points[hull.vertices,0])
    tempy = np.mean(hull.points[hull.vertices,1])
    tempz = np.mean(hull.points[hull.vertices,2])
    x = hull.points[hull.vertices,0]
    y = hull.points[hull.vertices,1]
    z = hull.points[hull.vertices,2] 
    hullx.append(hull.points[hull.vertices,0])
    hully.append(hull.points[hull.vertices,1])
    hullz.append(hull.points[hull.vertices,2])   
    # print ("Number of convex hull vertices",len(x),len(y),len(z))
    # print ("Convex hull points", x,y,z)  
    
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # ax.set_xlabel('x1')
    # ax.set_ylabel('x2')
    # ax.set_zlabel('x3')
    # ax.scatter3D(x, y, z, zdir='z',cmap='viridis')
    # ax.scatter3D(cx, cy, cz)
    # plt.show()
    
    #attempt
    # P = hull.points
    # x =[float(item[0]) for item in P]
    # y =[float(item[1]) for item in P]
    # z = [float(item[2]) for item in P]
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.2)
    # #plt.show()

    # print ("Centroid points", cx,cy,check_angle(cz))  
    convex_points = []
    convex_points = [list(a) for a in zip(x,y,z)]
    # print (convex_points)
    #input("Press Enter to continue...")
    #return (tempx,tempy,check_angle(tempz),convex_points)
    return (tempx,tempy,tempz,convex_points)

def maxdistance_set_target(est_set):
    max_distance = 0 
    for vertices in est_set:
        dist = math.sqrt((vertices[0]-x1t)**2 + (vertices[1]-x2t)**2)
        if dist > max_distance: max_distance = dist
    return max_distance

def translational_convergence(max_distance, x_under,x_over,y_under,y_over):
    offloading_condition = 0 
    vertices = list(itertools.product([x_under,x_over],[y_under,y_over]))
    for point in vertices:
        dist = math.sqrt((point[0]-x1t)**2 + (point[1]-x2t)**2)
        if dist >= max_distance: 
            print ("oleeee OFFLOAD")
            offloading_condition = 1
    return offloading_condition

def network_overhead_calc(x1,x2):
    min_dist = 10000
    for vertices in ap_locations:
        dist = math.sqrt((vertices[0]-x1)**2 + (vertices[1]-x2)**2)
        if dist < min_dist: 
            min_dist = dist  
    print ("Distance from AP", dist)
    c_bytes = network_delay_distance[int(min_dist)]
    print ("c bytes", c_bytes)
    # random uniform for image megabytes
    image_bytes = random.uniform(0.05,0.1) * 1000000
    network_overhead = image_bytes / c_bytes
    print ("network overhead", network_overhead)
    #input("Press Enter to continue...")
    return network_overhead

def computing_overhead_calc():
    mu, sigma = 0.75, 0.16 # mean and standard deviation
    cores = np.random.normal(mu, sigma)
    if cores < 0.25: cores = 0.25
    if cores > 1.5: cores =  1.5
    cores= random.uniform(0.25,2)
    computing_overhead = -1.34*cores+3.675
    print ("computing overhead", computing_overhead)
    #input("Press Enter to continue...")
    return computing_overhead

def utility_function(offloading_condition,network_overhead,computing_overhead,max_distance):
    if offloading_condition == 1: return 20 
    convergence = estimation_volume[-1] / max_distance
    c1 = 12
    c2 = 2
    utility  = c1 * convergence - c2 * (network_overhead+computing_overhead) + 10
    print ("vol/d:",  convergence)
    print ("e_uoff:",  (network_overhead+computing_overhead))
    print ("utility", utility)
    filename = "./utility.csv"
    with open(filename, 'a') as myfile:
	    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
	    # If opened for the first time, insert header row           
	    if os.path.getsize(filename) == 0:
	        wr.writerow(["max_distance_from_target","volume_estimation_set","vol/dmax","network_overhead","computing_overhead","e_offl","utility"])
	    wr.writerow([max_distance,estimation_volume[-1],convergence,network_overhead,computing_overhead,-(network_overhead+computing_overhead),utility])
    #input("Press Enter to continue...")
    return utility 

def simulation():
    starting_time = time.time()
    number_offloading = 0 
    #initialization of simulation
    if TRAN == True:
        x1in = x1 = 5#random.randint(0,25)
        x2in = x2 = 5#random.randint(0,25) 
        x3in = x3 = 1#random.randrange(100,128)/100 #random.randrange(10,628)/100
    else : 
        x1in = x1 = random.randint(0,15)
        while (x1 == x1t):
            x1in = x1 = random.randint(0,15)
        x2in = x2 = random.randint(0,15)
        while (x2 == x2t):
            x2in = x2 = random.randint(0,15) 
        x3in = x3 = random.randrange(100,128)/100 #random.randrange(10,628)/100
    xmin = x1 - d1_camera
    xmax = x1 + d1_camera
    ymin = x2 - d2_camera
    ymax = x2 + d2_camera
    x3min = x3 - d3_camera
    x3max = x3 + d3_camera 
    est_set = [[x1-d1_camera,x1+d1_camera],[x2-d2_camera,x2+d2_camera],[x3-d3_camera,x3+d3_camera]]
    xvalue,yvalue,x3value,cube = calculateEstimationSet(est_set[0][0],est_set[0][1], est_set[1][0], est_set[1][1], est_set[2][0],est_set[2][1])
    est_set = cube
    estimation_parallelotope.append(list(itertools.product([x1-d1_camera,x1+d1_camera],[x2-d2_camera,x2+d2_camera],[x3-d3_camera,x3+d3_camera])))
    initial_set = ["["+str(x1- d1_camera)+","+str(x1+d1_camera)+"]","["+str(x2- d2_camera)+","+str(x2+d2_camera)+"]","["+str(x3- d3_camera)+","+str(x3+d3_camera)+"]",]
    flow_volume.append(0)
    estimation_volume.append(0)
    steps = 0 
    steps2 = 0 
    changes_in_m = 0 
    motion = 0              # 0 for translational , 1 for rotational 
    distance = math.sqrt((x1-x1t)**2 + (x2-x2t)**2)
    print ("Initial Position:", x1,x2,"{:10.0f}".format(math.degrees(x3)))
    print ("Target Position: ", x1t,x2t)
    print ("\n")
    offload = True
    while (distance>e):
        #input("Press Enter to continue...")
        df = calculateangle (x1,x2,x3) 
        #TODO fix zero division error
        if abs(x1-x1t)/(x1-x1t) >= 0 or abs(x2-x2t)/(x2-x2t) >= 0 :
            angle = x3 - math.pi - math.atan2((x2-x2t),(x1-x1t)) 
        else:
            angle = x3 - abs(math.atan((x2-x2t)/(x1-x1t)))  
        #m = abs(angle) * (math.sqrt((x1-x1t)**2 + (x2-x2t)**2))
        m = abs(df) * (math.sqrt((x1-x1t)**2 + (x2-x2t)**2))   # check if we can proceed to rotational motion
        #TODO implement new M set
        distance = (math.sqrt((x1-x1t)**2 + (x2-x2t)**2))
        d = d3_rot[0] # implement the same for [1]
        l = math.sqrt(2/d**2*abs(1-math.cos(d)))
        u1 = 4*distance/l * math.cos(d+angle) 
        #print ("u are:", distance,l,(d+angle), math.cos(d+angle))
        d = d3_rot[1] # implement the same for [1]
        l = math.sqrt(2/d**2*abs(1-math.cos(d)))
        u2 = 4*distance/l * math.cos(d+angle) 
        #print ("u are:", distance,l,(d+angle), math.cos(d+angle))
        u = min(u1,u2)
        print (u)
        #print ("u are:", u1,u2)
        #input("Press Enter to continue...")
        #m = run_time**2*u * math.sin(2*x3+2*d3_rot[1])
        #m2 = 4*math.pow(d,4)*(math.cos(d3_tran[1]+angle))**2
       
        if u > 0 and offload :
            if motion != 0 : 
                motion = 0
                changes_in_m +=1 
            
            #x1,x2,x3,gainx1,gainx2,gainx3,distance = translational(x1,x2,x3)
            #u = translational_noisy(est_set)
            gainx1 = u *discretization
            gainx2 = u *discretization
            gainx3 = 0  
            x1f = x1 +  gainx1 * math.cos(x3) 
            x2f = x2 +  gainx2 * math.sin(x3) 
            distance_traversed = (x1f-x1)**2 + (x2f-x2)**2

            d1 = [d1_tran[0], d1_tran[1]]
            d2 = [d2_tran[0], d2_tran[1]]
            d3 = [d3_tran[0] ,d3_tran[1]]
            mode = "tran"
            jumptime = 1 
        else:
            if offload == False: offload = True
            if motion != 1 : 
                motion = 1
                changes_in_m +=1 
            u, x3 = rotational_noisy(est_set)
            gainx3 = u #* r/l * k1
            gainx1 = 0 
            gainx2 = 0 
            #x1,x2,x3,gainx1,gainx2,gainx3 = rotational(x1,x2,x3,df)
            d1 = d1_rot
            d2 = d2_rot 
            d3 = [d3_rot[0], d3_rot[1]]
            distance_traversed = 0 

            mode = "rot"
            jumptime = 0 
        
        x_under,x_over,y_under,y_over,x3_under,x3_over = run_flow(time,d1,d2,d3,gainx1,gainx2,gainx3,jumptime,initial_set,mode)
        max_distance = maxdistance_set_target(est_set)
        offloading_condition = 0 
        if mode=="tran":
            offloading_condition = translational_convergence(max_distance,x_under,x_over,y_under,y_over)
        
        network_overhead = network_overhead_calc(x1,x2)
        computing_overhead = computing_overhead_calc()
        remote_time = int((network_overhead+computing_overhead)/0.1)
        utility = utility_function(offloading_condition,network_overhead,computing_overhead,max_distance)

        if (steps >break_steps):
            break
        if (utility > 6.5 or steps2 > camera_steps):    # use camera and make estimation set a point! 
            #input("Press Enter to continue...")
            offload = False
            number_offloading += 1
            #random point in estimation set! not just the representative
            listx = [item[0] for item in est_set]
            listy = [item[1] for item in est_set]
            listx3 = [item[2] for item in est_set]
            x1 = random.uniform(min(listx),max(listx))
            x2 = random.uniform(min(listy),max(listy))
            x3 = random.uniform(min(listx3),max(listx3))
            initial_set = ["["+str(x1- d1_camera)+","+str(x1+d1_camera)+"]","["+str(x2- d2_camera)+","+str(x2+d2_camera)+"]","["+str(x3- d3_camera)+","+str(x3+d3_camera)+"]",]
            xmin = x1 - d1_camera
            xmax = x1 + d1_camera
            ymin = x2 - d2_camera
            ymax = x2 + d2_camera
            x3min = x3 - d3_camera
            x3max = x3+ d3_camera 
            max_distance = maxdistance_set_target(est_set)
            x1,x2,x3,est_set = calculateEstimationSet(xmin,xmax,ymin,ymax,x3min,x3max)
            initial_set, xmin,xmax,ymin,ymax,x3min,x3max = get_initialset(est_set)

            steps2 = 0 
            steps+=remote_time # simulate transmission + plus computation time! unfortunately this is in sec
            for i in range (remote_time):
                lists_renew(x1,x2,x3,max_distance,m,df,angle)
                flow_volume.append(flow_volume[-1])
                estimation_volume.append(estimation_volume[-1])
        else:
            x1,x2,x3,est_set =compatibility(time,d1,d2,d3,gainx1,gainx2,gainx3,jumptime,initial_set,mode,xmin,xmax,ymin,ymax,x3min,x3max,distance_traversed,x3,x_under,x_over,y_under,y_over,x3_under,x3_over)
            initial_set, xmin,xmax,ymin,ymax,x3min,x3max = get_initialset(est_set)
            lists_renew(x1,x2,x3,max_distance,m,df,angle) 
            steps +=1 
            steps2 += 1 
        # collision , closest_point = check_feasibility(est_set)
        # if collision:
        #     x1 = closest_point[0]
        #     x2 = closest_point[1]
        #     x3 = closest_point[2]

        #TODO if we want to assign the representative to a critical vertice (close to obstacle this is the place)
        #x1,x2,x3,est_set = calculateConvexHull(est_set,vertices)

        #input("Press Enter to continue...")
        print ("Final Position: ","{:10.2f}".format(x1t),"{:10.2f}".format(x2t), x3)
        print ("Representative Position: ","{:10.2f}".format(x1),"{:10.2f}".format(x2), "{:10.2f}".format(x3))
        print ("M: ",m)

        #distance = math.sqrt((x1-x1t)**2 + (x2-x2t)**2)
        print ("Distance difference:", "{:10.2f}".format(distance))
        print ("Angle difference:", "{:10.2f}".format(df))
        print ("Steps: ",steps)
        print ("Volume of estimation set: ", estimation_volume[-1])
        print ("number_offloading ", number_offloading)
        print ("\n ")
    #end of while 


    print ("Initial Position:", x1in,x2in,"{:10.0f}".format(math.degrees(x3in)))
    print ("Target Position: ", x1t,x2t)
    print ("Final Position: ","{:10.2f}".format(x1),"{:10.2f}".format(x2),"{:10.0f}".format(math.degrees(x3)))
    print ("Distance error tolerance: ", e)
    print ("M tolerance: ", M)
    print ("Number of changes in M:", changes_in_m)
    print ("Steps:", steps)
    print ("Real time needed for simulation in sec :", "{:10.2f}".format(time.time() - starting_time))
    plots(steps)  
    #plot_cube()
    
if __name__ == '__main__':
    simulation()
