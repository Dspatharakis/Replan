#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import time 
import subprocess
import itertools
import cvxpy as cp
from pypoman import compute_polytope_halfspaces
from numpy import exp,arange , arctan , sqrt
from scipy import signal
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from subprocess import call
from matplotlib import cm

# some "global" variables
plot = True     # True for plotting
r = 6.6 /2.0    # 6.6(cm) radius of wheels in cm 
l = 13.2        # (cm) distance between the two weels 
x1t = random.randint(0,25)  # starting position for x1
x2t = random.randint(0,25)  # starting position for x2
M = 5           # c boundary for M set
e = 1           # break condition for distance from target position
k1 = 10.8#10.8  # control gain for rotational motion
k2 = 0.7#0.05   # control gain for x1 in translational motion
k3 = 0.7#0.05   # control gain for x2 in translational motion

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

# some variables
run_time = 0.01      # discretization time 
break_steps = 400    # break time! 
plotvar = 25         # define after how many steps the cube will be plotted! 
camera_steps = 150   # define after how many steps we use the camera. Later on this maybe regarding the volume of the box
# static disturbances
d1 = [-0.01,0.01]
d2 = [-0.01,0.01]
d3 = [-0.1,0.1]



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
        plt.show()

    # x1
    plt.plot(t, x1_set, 'g.',label='x1 over time')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('x1 (cm)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("./output/x1.png")
        plt.show()

    # x2
    plt.plot(t, x2_set, 'b.',label='x2 over time')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('x2 (cm)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("./output/x2.png")
        plt.show()

    # x3
    plt.plot(t, x3_set, 'r.',label='x3 over time')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('x3 (degrees)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("./output/x3.png")
        plt.show()

    # m
    plt.plot(t, m_set, 'b.',label='m over time/ limit=%d' %M)
    plt.axhline(M,color='red')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('m ', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("./output/m.png")
        plt.show()

    # d
    plt.plot(t, d_set, 'b.',label='d over time/ limit=%d' %e)
    plt.axhline(e,color='red')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('d (cm)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("./output/d.png")
        plt.show()

    # angle dif 
    plt.plot(t, angledf_set, 'r.',label='angledif over time')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('angledif  (degrees)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("./output/angledif.png")
        plt.show()

    # df
    plt.plot(t, df_set, 'r.',label='angle needed to rotate in order to look target')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('angle needed to rotate in order to look target (degrees)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("./output/df.png")
        plt.show()

    # difference between two methods
    plt.plot(t, temp_set, 'r.',label='difference between two methods')
    plt.xlabel('steps', color='#1C2833')
    plt.ylabel('difference between two methods (degrees)', color='#1C2833')
    plt.legend(loc='upper left')
    plt.grid()
    if plot:
        plt.savefig("./output/methods.png")
        plt.show()

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
        printx = []
        printy = []
        printx3 = []
        printx.append(hullx[o][0])      
        printx.append(hullx[o][1])
        printy.append(hully[o][0])
        printy.append(hully[o][1])
        printx3.append(hullz[o][0])
        printx3.append(hullz[o][1])
        #sthn plot pane ta x values duo shmeiwn kai meta ta y values 2 shmeiwn
        ax.plot(printx, printy, printx3, color=   'g')
        # 
        printx = []
        printy = []
        printx3 = []
        printx.append(hullx[o][2])      
        printx.append(hullx[o][3])
        printy.append(hully[o][2])
        printy.append(hully[o][3])
        printx3.append(hullz[o][2])
        printx3.append(hullz[o][3])
        #sthn plot pane ta x values duo shmeiwn kai meta ta y values 2 shmeiwn
        ax.plot(printx, printy, printx3, color=   'g')
        # 
        printx = []
        printy = []
        printx3 = []
        printx.append(hullx[o][4])      
        printx.append(hullx[o][5])
        printy.append(hully[o][4])
        printy.append(hully[o][5])
        printx3.append(hullz[o][4])
        printx3.append(hullz[o][5])
        #sthn plot pane ta x values duo shmeiwn kai meta ta y values 2 shmeiwn
        ax.plot(printx, printy, printx3, color=   'g')
        # 
        printx = []
        printy = []
        printx3 = []
        printx.append(hullx[o][6])      
        printx.append(hullx[o][7])
        printy.append(hully[o][6])
        printy.append(hully[o][7])
        printx3.append(hullz[o][6])
        printx3.append(hullz[o][7])
        #sthn plot pane ta x values duo shmeiwn kai meta ta y values 2 shmeiwn
        ax.plot(printx, printy, printx3, color=   'g')
        # 
        printx = []
        printy = []
        printx3 = []
        printx.append(hullx[o][0])      
        printx.append(hullx[o][2])
        printy.append(hully[o][0])
        printy.append(hully[o][2])
        printx3.append(hullz[o][0])
        printx3.append(hullz[o][2])
        #sthn plot pane ta x values duo shmeiwn kai meta ta y values 2 shmeiwn
        ax.plot(printx, printy, printx3, color=   'g')
        # 
        printx = []
        printy = []
        printx3 = []
        printx.append(hullx[o][1])      
        printx.append(hullx[o][3])
        printy.append(hully[o][1])
        printy.append(hully[o][3])
        printx3.append(hullz[o][1])
        printx3.append(hullz[o][3])
        #sthn plot pane ta x values duo shmeiwn kai meta ta y values 2 shmeiwn
        ax.plot(printx, printy, printx3, color=   'g')
        # 
        printx = []
        printy = []
        printx3 = []
        printx.append(hullx[o][0])      
        printx.append(hullx[o][4])
        printy.append(hully[o][0])
        printy.append(hully[o][4])
        printx3.append(hullz[o][0])
        printx3.append(hullz[o][4])
        #sthn plot pane ta x values duo shmeiwn kai meta ta y values 2 shmeiwn
        ax.plot(printx, printy, printx3, color=   'g')
        # 
        printx = []
        printy = []
        printx3 = []
        printx.append(hullx[o][1])      
        printx.append(hullx[o][5])
        printy.append(hully[o][1])
        printy.append(hully[o][5])
        printx3.append(hullz[o][1])
        printx3.append(hullz[o][5])
        #sthn plot pane ta x values duo shmeiwn kai meta ta y values 2 shmeiwn
        ax.plot(printx, printy, printx3, color=   'g')
        # 
        printx = []
        printy = []
        printx3 = []
        printx.append(hullx[o][4])      
        printx.append(hullx[o][6])
        printy.append(hully[o][4])
        printy.append(hully[o][6])
        printx3.append(hullz[o][4])
        printx3.append(hullz[o][6])
        #sthn plot pane ta x values duo shmeiwn kai meta ta y values 2 shmeiwn
        ax.plot(printx, printy, printx3, color=   'g')
        # 
        printx = []
        printy = []
        printx3 = []
        printx.append(hullx[o][5])      
        printx.append(hullx[o][7])
        printy.append(hully[o][5])
        printy.append(hully[o][7])
        printx3.append(hullz[o][5])
        printx3.append(hullz[o][7])
        #sthn plot pane ta x values duo shmeiwn kai meta ta y values 2 shmeiwn
        ax.plot(printx, printy, printx3, color=   'g')
        # 
        printx = []
        printy = []
        printx3 = []
        printx.append(hullx[o][2])      
        printx.append(hullx[o][6])
        printy.append(hully[o][2])
        printy.append(hully[o][6])
        printx3.append(hullz[o][2])
        printx3.append(hullz[o][6])
        #sthn plot pane ta x values duo shmeiwn kai meta ta y values 2 shmeiwn
        ax.plot(printx, printy, printx3, color=   'g')
        # 
        printx = []
        printy = []
        printx3 = []
        printx.append(hullx[o][3])      
        printx.append(hullx[o][7])
        printy.append(hully[o][3])
        printy.append(hully[o][7])
        printx3.append(hullz[o][3])
        printx3.append(hullz[o][7])
        #sthn plot pane ta x values duo shmeiwn kai meta ta y values 2 shmeiwn
        ax.plot(printx, printy, printx3, color=   'g')



    #ax.scatter(printx, printy, printx3, c='r',s=100)
    #ax.plot(printx, printy, printx3, color='r')

    #ax.scatter3D(printx, printy, printx3, zdir='z',cmap='viridis')
    cxlist = []
    cylist = []
    cx3list = []
    for x in range (0 , len(hullx),plotvar): 
        cxlist.append(cx[x])      # ta x ana 10 kai ana shmeio tou kubou 
        cylist.append(cy[x])
        cx3list.append(cz[x])
    
    ax.scatter3D(cxlist, cylist, cx3list)
    #ax.scatter3D(cx, cy, cz)
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

def build_model(time, d1, d2, d3, u1, u2, u3, jumptime, initial_set, mode,output,name):
    with open('replan.model', 'w') as the_file:
        the_file.write('hybrid reachability\n''{\n''  state var x,y,x3,t \n'' setting\n''{\n'
        'fixed steps '+str(run_time)+'\n''time '+str(run_time)+'\n' # time  
        'remainder estimation 1e-2\n''identity precondition\n''gnuplot octagon '+ output+'\n'
        'fixed orders 8\n''cutoff 1e-15\n''precision 2000\n''output '+name+'\n''max jumps 1\n''print on\n''}\n')
        the_file.write('modes\n'
        '{\n''tran\n''{ \n''nonpoly ode\n''{\n'"t' = 1\n")
        the_file.write("x3' ="+ str(d3) + "\n" # to d 1
        "x' = "+str(u1)+"*cos(x3) \n"   # to u 1
        "y' = "+str(u2)+"*sin(x3) \n""}\n") # to u 2
        the_file.write("inv\n""{\n""t <= 1""\n""t<=1""\n") 
        # for item in invtran:
        #     the_file.write(item+"\n")
        the_file.write("\n"   
        "}\n""}\n"
        "rot\n""{\n""nonpoly ode\n""{\n"
        "x3' = "+str(u3)+"\n"  # to u3
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
    return octagon_vertices[:-1]

def get_initialset(est_set):
    listx = [item[0] for item in est_set]
    listy = [item[1] for item in est_set]
    listx3 = [item[2] for item in est_set]
    initial_set = ["["+str(min(listx))+","+str(max(listx))+"]","["+str(min(listy))+","+str(max(listy))+"]","["+str(min(listx3))+","+str(max(listx3))+"]",]
    # print ("Initial set: ")
    # print (initial_set)
    # input("Press Enter to continue...")
    return initial_set, min(listx),max(listx), min(listy),max(listy), min(listx3),max(listx3)

def run_flow(time, d1, d2, d3, u1, u2, u3, jumptime , initial_set, mode,xmin,xmax,ymin,ymax,x3min,x3max,distance,x3):
    output = "x,t"
    build_model(time, d1, d2, d3, u1/run_time, u2/run_time, u3/run_time, jumptime , initial_set, mode , output, "x") 
    rc = call(["./run.sh","x"])
    verticesx = get_numbers()
    #input("Press Enter to continue...")
    output = "y,t"
    build_model(time, d1, d2, d3, u1/run_time, u2/run_time, u3/run_time,  jumptime, initial_set, mode , output, "y") 
    rc = call(["./run.sh","y"])
    verticesy = get_numbers()
    #input("Press Enter to continue...")
    output = "x3,t"
    build_model(time, d1, d2, d3, u1/run_time, u2/run_time, u3/run_time,jumptime, initial_set, mode , output, "x3") 
    rc = call(["./run.sh","x3"])
    verticesx3 = get_numbers()
    print (mode)
    print ("Gains: ",u1,u2,u3)
    #input("Press Enter to continue...")
    listx = [item[0] for item in verticesx]
    listy = [item[0] for item in verticesy]
    listx3 = [item[0] for item in verticesx3]
    x_under = min(listx)
    x_over = max(listx)
    y_under = min(listy)
    y_over = max(listy)
    x3_under = min(listx3)
    x3_over = max(listx3)
    print ("Flow:",x_under,x_over,y_under,y_over,x3_under,x3_over)
    measur1 = distance
    w =  2.2 #random.uniform(0.001,0.002)
    measur2 = u3
    w2 = 0.1* abs(measur2) #random.uniform(0.001,0.002)
    yrootmin, yrootmax = optimize_quadratic(y_under,y_over, ymin,ymax,measur1,w)
    xrootmin, xrootmax = optimize_quadratic(x_under,x_over, xmin,xmax,measur1,w)

    xminstar = max(x_under, (xmin))
    xmaxstar = min(x_over, (xmax +math.sqrt(measur1+w)))
    yminstar = max(y_under, (ymin))
    ymaxstar = min(y_over, (ymax +math.sqrt(measur1+w)))
    #print ("Measur:" , measur1,measur2)
    x3minstar, x3maxstar = optimization_problem(x3_under,x3_over,x3min,x3max,measur2,w2)    
    print ("Optimization", xminstar,xmaxstar,yminstar,ymaxstar,x3minstar,x3maxstar)
    #vertices = list(itertools.product([xminstar,xmaxstar],[yminstar,ymaxstar],[x3minstar,x3maxstar]))
    xvalue,yvalue,x3value,cube = calculateEstimationSet(xminstar,xmaxstar,yminstar,ymaxstar,x3minstar,x3max )
    #return vertices 
    return xvalue,yvalue,x3value,cube

def optimization_problem(v3min,v3max,z3min,z3max,y,w):
    #TODO may need x_3 to be in [0,6.7] instead ofhaving negative values. CHECK!
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
    v = cp.Variable()
    z = cp.Variable()
    prob1 = cp.Problem(cp.Minimize(cp.square(v-z)),
                        [cp.square(v-z)   <= y+w , 
                        z >= zmin , z <= zmax, 
                        v >= vmin , v <= vmax ])
    
    try:
        minimum = prob1.solve()
    except Exception as e:
        print(e)

    v = cp.Variable()
    z = cp.Variable()
    prob1 = cp.Problem(cp.Maximize(v-z),
                        [cp.square(v-z)   <= y-w , 
                        z >= zmin , z <= zmax, 
                        v >= vmin , v <= vmax ])
    
    try:
        maximum = prob1.solve()
    except Exception as e:
        print(e)

    return round(minimum,8) , round(maximum,8)

def calculateEstimationSet(xmin,xmax,ymin,ymax,x3min,x3max):  
    #TODO append cx,cy,cz and hullx,hully,hullz
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

def simulation():
    #initialization of simulation
    x1in = x1 = random.randint(0,25)
    x2in = x2 = random.randint(0,25) 
    x3in = x3 = random.randrange(100,128)/100 #random.randrange(10,628)/100
    xmin = xmax = x1
    ymin = ymax = x2
    x3min = x3max = x3
    est_set = [[x1,x2,x3]]
    initial_set = ["["+str(x1)+","+str(x1)+"]","["+str(x2)+","+str(x2)+"]","["+str(x3)+","+str(x3)+"]",]
    steps = 0 
    steps2 = 0 
    changes_in_m = 0 
    motion = 0              # 0 for translational , 1 for rotational 
    d = math.sqrt((x1-x1t)**2 + (x2-x2t)**2)
    print ("Initial Position:", x1,x2,"{:10.0f}".format(math.degrees(x3)))
    print ("Target Position: ", x1t,x2t)
    print ("\n")
    while (d>e):
        #input("Press Enter to continue...")
        df = calculateangle (x1,x2,x3) 
        #TODO fix zero division error
        if abs(x1-x1t)/(x1-x1t) >= 0 or abs(x2-x2t)/(x2-x2t) >= 0 :
            angle = x3 - math.pi - math.atan2((x2-x2t),(x1-x1t)) 
        else:
            angle = x3 - abs(math.atan((x2-x2t)/(x1-x1t)))  
        #m = abs(angle) * (math.sqrt((x1-x1t)**2 + (x2-x2t)**2))
        m = abs(df) * (math.sqrt((x1-x1t)**2 + (x2-x2t)**2))   # check if we can proceed to rotational motion
        if m < M :
            if motion != 0 : 
                motion = 0
                changes_in_m +=1 
            x1,x2,x3,gainx1,gainx2,gainx3,distance = translational(x1,x2,x3)
            mode = "tran"
            jumptime = 1 
        else:
            if motion != 1 : 
                motion = 1
                changes_in_m +=1 
            x1,x2,x3,gainx1,gainx2,gainx3 = rotational(x1,x2,x3,df)
            distance = 0 
            mode = "rot"
            jumptime = 0 
        x1,x2,x3,est_set = run_flow(time,d1,d2,d3,gainx1,gainx2,gainx3,jumptime,initial_set,mode,xmin,xmax,ymin,ymax,x3min,x3max,distance,x3)
        #x1,x2,x3,est_set = calculateConvexHull(est_set,vertices)
        initial_set, xmin,xmax,ymin,ymax,x3min,x3max = get_initialset(est_set)
        print ("Final Position: ","{:10.2f}".format(x1t),"{:10.2f}".format(x2t), x3)
        print ("Representative Position: ","{:10.2f}".format(x1),"{:10.2f}".format(x2), x3)
        print ("M: ",m)
        d = math.sqrt((x1-x1t)**2 + (x2-x2t)**2)
        print ("Distance difference:", d)
        print ("Angle difference:", df)
        print ("Steps: ",steps)
        print ("\n ")
        lists_renew(x1,x2,x3,d,m,df,angle) 
        steps +=1 
        steps2 += 1 
        if steps > break_steps : break
        if steps2 >= camera_steps :    # use camera and make estimation set a point! 
            initial_set = ["["+str(x1)+","+str(x1)+"]","["+str(x2)+","+str(x2)+"]","["+str(x3)+","+str(x3)+"]",]
            xmin = xmax = x1
            ymin=ymax = x2
            x3min = x3max = x3
            steps2 = 0 
    print ("Initial Position:", x1in,x2in,"{:10.0f}".format(math.degrees(x3in)))
    print ("Target Position: ", x1t,x2t)
    print ("Final Position: ","{:10.2f}".format(x1),"{:10.2f}".format(x2),"{:10.0f}".format(math.degrees(x3)))
    print ("Distance error tolerance: ",e)
    print ("M tolerance: ",M)
    print ("Number of changes in M:",changes_in_m)
    print ("Steps:",steps)
    plots(steps)  
    plot_cube()
    
if __name__ == '__main__':
    simulation()
