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

plot = True # True for plotting
r = 6.6 /2.0 # 6.6(cm) radius of wheels in cm 
l = 13.2     # (cm) distance between the two weels 
T = 0.3      #discretization time 
x1t = random.randint(0,250)  # starting position for x1
x2t = random.randint(0,250)  # starting position for x2
M = 5    # c boundary for M set
e = 5    # break condition for distance from target position
k1 = 10.8 # control gain for rotational motion
k2 = 0.1 # control gain for x1 in translational motion
k3 = 0.1 # control gain for x2 in translational motion
# some initialization of lists for printing
x1_set = []
x2_set = []
x3_set = [] 
d_set = []
m_set = []
df_set = []
angledf_set = []
temp_set =[]
Χ = []

def rotational(x1,x2,x3,df):
    gain = T * r/l * k1 * df  # (x3 - math.atan2((x2-x2t),(x1-x1t))) 
    #print (gain)
    x1 = x1 #+  100 * random.uniform(-gain,gain)   #TODO 100 is very random right now
    x2 = x2 #+  100 * random.uniform(-gain,gain) 
    x3 = x3 + gain #+ 5 * random.uniform(-gain,gain) 
    return x1,x2,x3,0,0,gain
    
def translational(x1,x2,x3):
    gainx1 = T * r * k2 * math.cos(x3) * math.sqrt((x1-x1t)**2+(x2-x2t)**2)
    gainx2 = T * r * k3 * math.sin(x3) * math.sqrt((x1-x1t)**2+(x2-x2t)**2)
    #print ("gainx1 : ", gainx1,gainx2)
    x1 = x1 +  gainx1 #+ 0.4 * random.uniform(-gainx1,gainx1) 
    x2 = x2 +  gainx2 #+ 0.4 * random.uniform(-gainx2,gainx2)
    x3 = x3 #+  0.01 * random.uniform(-math.sqrt(gainx1**2+gainx2**2),math.sqrt(gainx1**2+gainx2**2))
    return x1,x2,x3, gainx1, gainx2, 0

def check_angle(x3):
    if x3 > 2*math.pi : 
        x3 = x3 - 2*math.pi
    if x3 < 0:
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
        plt.savefig("grid.png")
        plt.show()

    # x1
    plt.plot(t, x1_set, 'g.',label='x1 over time')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('x1 (cm)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("x1.png")
        plt.show()

    # x2
    plt.plot(t, x2_set, 'b.',label='x2 over time')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('x2 (cm)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("x2.png")
        plt.show()

    # x3
    plt.plot(t, x3_set, 'r.',label='x3 over time')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('x3 (degrees)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("x3.png")
        plt.show()

    # m
    plt.plot(t, m_set, 'b.',label='m over time/ limit=%d' %M)
    plt.axhline(M,color='red')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('m ', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("m.png")
        plt.show()

    # d
    plt.plot(t, d_set, 'b.',label='d over time/ limit=%d' %e)
    plt.axhline(e,color='red')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('d (cm)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("d.png")
        plt.show()

    # angle dif 
    plt.plot(t, angledf_set, 'r.',label='angledif over time')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('angledif  (degrees)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("angledif.png")
        plt.show()

    # df
    plt.plot(t, df_set, 'r.',label='angle needed to rotate in order to look target')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('angle needed to rotate in order to look target (degrees)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("df.png")
        plt.show()

    # difference between two methods
    plt.plot(t, temp_set, 'r.',label='difference between two methods')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('difference between two methods (degrees)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("methods.png")
        plt.show()

def lists_renew(x1,x2,x3,d,m,df,angle):
    x1_set.append(x1)
    x2_set.append(x2)
    x3_set.append(math.degrees(x3))
    d_set.append(d)
    m_set.append(m)
    df_set.append(math.degrees(df))
    angledf_set.append(math.degrees(angle))
    temp_set.append(math.degrees(df)-math.degrees(angle))

def propaggation_of_set(est_set,gainx1,gainx2,gainx3,y1,y2,y3):  
    temp = []
    temp_set = set()
    start_time = time.time()
    non_compatible = 0 
    print ("Gains: ",gainx1,gainx2,gainx3)
    #measurement noise 
    w1 = random.uniform(0.1,2.0)
    w2 = random.uniform(0.1,2.0)
    w3 = random.uniform(0.1,2.0)
    w4 = random.uniform(0.1,2.0)
    w5 = random.uniform(0.1,0.2)
    w6 = random.uniform(0.1,0.2)

    #process noise 
    o1 = random.uniform(0.1,1.0)
    o2 = random.uniform(0.1,1.0)
    o3 = random.uniform(0.1,1.0)
    o4 = random.uniform(0.1,1.0)
    o5 = random.uniform(0.01,0.2)
    o6 = random.uniform(0.01,0.2)
    #print ("os: ",o1,o2,o3,o4,o5,o6)

    for item in est_set:
        #print (item)
        x1=item[0]
        x2=item[1]
        x3=item[2]
        # for d1 in np.arange(-1,1.9,0.5):   
        #     for d2 in np.arange(-1,1.9,0.5):
        #         for d3 in np.arange(-0.1,0.19,0.1):
        for d1 in np.arange(-o1,o2,0.1):   
            for d2 in np.arange(-o3,o4,0.1):
                for d3 in np.arange(-o5,o6,0.1):
                    x33 = check_angle(x3+gainx3+d3)
                    x33 = round(x33,2)
                    x11 = round(x1+gainx1+d1,1)
                    x22 = round(x2+gainx2+d2,1)
                    x11 = check_grid(x11)
                    x22 = check_grid(x22)
                    #print (x33,check_angle(y3))
                    #print (x11,y1)
                    #print (x22,y2)
                    if ((x11,x22,x33) not in temp_set) :
                        if ( (x11-x1-gainx1 < w1) and (-x11+x1+gainx1 < w2 ) and (x22-x2-gainx2 < w3) and (-x22+x2+gainx2 < w4) and (x33-check_angle(x3+gainx3) < w5) and (-x33+check_angle(x3+gainx3) < w6) ):
                            temp.append([x11,x22,x33])
                            temp_set.add((x11,x22,x33))
                        else: non_compatible += 1 

    #print (temp)
    print ("Not compatible ", non_compatible) 
    #print ("w:", w1,w2,w3,w4,w5,w6)
    print ("Number of possible positions :", len(temp))
    print ("Time for iterations", time.time()-start_time)
    start_time = time.time()
    points = np.vstack(temp)
    print ("Time for vstack", time.time()-start_time)
    start_time = time.time()
    hull = ConvexHull(points)#,qhull_options='QJ')
    print ("Number of points :", len(points))
    print ("Time for convex hull", time.time()-start_time)
    #Get centoid
    cx = np.mean(hull.points[hull.vertices,0])
    cy = np.mean(hull.points[hull.vertices,1])
    cz = np.mean(hull.points[hull.vertices,2])
    x = hull.points[hull.vertices,0]
    y = hull.points[hull.vertices,1]
    z = hull.points[hull.vertices,2]   
    print ("Number of convex hull vertices",len(x),len(y),len(z))
    #print ("Convex hull points", x,y,z)  
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # ax.set_xlabel('x1')
    # ax.set_ylabel('x2')
    # ax.set_zlabel('x3')
    # ax.scatter3D(x, y, z, zdir='z',cmap='viridis')
    # ax.scatter3D(cx, cy, cz)
    # plt.show()
    print ("Centroid points", cx,cy,cz)  
    #input("Press Enter to continue...")
    convex_points = []
    convex_points = [list(a) for a in zip(x,y,z)]
    print (convex_points)
    # convex_points = overapproximation 
    return (cx,cy,cz,convex_points)
    #return (cx,cy,cz,temp)


def simulation():
    x1in = x1 = random.randint(0,250)
    x2in = x2 = random.randint(0,250) 
    x3in = x3 = random.randint(0,6)
    est_set = [[x1,x2,x3]]
    changes_in_m = 0 
    motion = 0 # 0 for translational , 1 for rotational 
    d = math.sqrt((x1-x1t)**2 + (x2-x2t)**2)
    print ("Initial Position:", x1,x2,"{:10.0f}".format(math.degrees(x3)))
    print ("Target Position: ", x1t,x2t)
    print ("\n")
    steps = 0 
    steps2 = 0 
    while (d>e):
        print ("Represantitive Positions: ",x1,x2,x3)
        df = calculateangle (x1,x2,x3) 
        ######## calculate m
        #TODO fix zero division error
        if abs(x1-x1t)/(x1-x1t) >= 0 or abs(x2-x2t)/(x2-x2t) >= 0 :
            angle = x3 - math.pi - math.atan2((x2-x2t),(x1-x1t)) 
        else:
            angle = x3 - abs(math.atan((x2-x2t)/(x1-x1t))) 
        ######## check if we can proceed to rotational motion
        m = abs(angle) * (math.sqrt((x1-x1t)**2 + (x2-x2t)**2))
        if m < M :
            if motion != 0 : 
                motion = 0
                changes_in_m +=1 
            x1,x2,x3,gainx1,gainx2,gainx3 = translational(x1,x2,x3)
        else:
            if motion != 1 : 
                motion = 1
                changes_in_m +=1 
            x1,x2,x3,gainx1,gainx2,gainx3 = rotational(x1,x2,x3,df)
        x1,x2,x3,est_set = propaggation_of_set(est_set,gainx1,gainx2,gainx3,x1,x2,x3)
        print ("Final Position: ","{:10.2f}".format(x1t),"{:10.2f}".format(x2t))
        print ("Steps: ",steps)
        print ("\n ")
        steps2 += 1 
        #if steps2 > 50:
        #    steps2 = 0
        #    est_set = [[x1,x2,x3]]
        x3 = check_angle(x3)
        d = math.sqrt((x1-x1t)**2 + (x2-x2t)**2)
        print ("Distance difference:", d)
        #print (x1-x1t)
        lists_renew(x1,x2,x3,d,m,df,angle) 
        steps +=1 
        if steps > 5000 : break

    print ("Initial Position:", x1in,x2in,"{:10.0f}".format(math.degrees(x3in)))
    print ("Target Position: ", x1t,x2t)
    print ("Final Position: ","{:10.2f}".format(x1),"{:10.2f}".format(x2),"{:10.0f}".format(math.degrees(x3)))
    print ("Number of changes in M:",changes_in_m)
    print ("Steps:",steps)
    plots(steps)


if __name__ == '__main__':
    simulation()

'''
#print ("d:", "{:10.4f}".format(d))
#print ("M: ", "{:10.4f}".format(m)) 
#print ("DF : ","{:10.4f}".format(df))
#print ("Angle difference: ","{:10.4f}".format( angle ))
#print ("{:10.4f}".format(x1),"{:10.4f}".format(x2),"{:10.4f}".format(x3))
#print ("Step",steps)
#print ("\n")

############################### plot 3d M 

def fun(x, y):
    return (x*y)
    #return ((np.arctan(y/x))*(np.sqrt(x**2+y**2)) ) 

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = np.arange(0.0, 25.0, 0.01)

y = np.arange(0.0, 25.0, 0.01)

z = np.arange(-60, 60, 0.01)
X, Y = np.meshgrid(x, y)

X,Y = np.where (X*Y <10 , X,0) , np.where (X*Y <10 , Y,0) 
X,Y = np.where (X*Y >-10 , X,0) , np.where (X*Y >-10 , Y,0) 

Z = fun(X,Y)

ax.plot_surface(X, Y, Z)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
plt.savefig('asd.png')
'''