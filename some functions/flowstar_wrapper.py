#!/usr/bin/env python3

#TODO input of function and parameterization of some lines 


a = "[-0.01,0.01]"
with open('replan.model', 'w') as the_file:
    the_file.write('hybrid reachability\n''{\n''  state var x,y,x3,t \n'' setting\n''{\n'
    'fixed steps 0.05\n''time 0.05\n' # time  
    'remainder estimation 1e-2\n''identity precondition\n''gnuplot octagon  x,y\n'
    'fixed orders 8\n''cutoff 1e-15\n''precision 2000\n''output ltv1_test\n''max jumps 1\n''print on\n''}\n')
    the_file.write('modes\n''{\n''tran\n''{ \n''nonpoly ode\n''{\n'"t' = 1\n")
    the_file.write("x3' ="+ a + "\n" # to d 
    "x' = 1*cos(x3) \n"   # to u 
    "y' = 1*sin(x3) \n""}\n") # to u 
    the_file.write("inv\n""{\n"
    "t >= 0\n"    # add compatibillity constraints 
    "}\n""}\n""rot\n""{\n""nonpoly ode\n""{\n"
    "x3' = 1\n"  # to u  
    "x' = [-0.1 , 0.1]\n""y' = [-0.1 , 0.1]\n" # to d
    "t' = 1\n""}\n""inv \n""{\n" # to d
 "x >= 0\n"  # add compatibillity constraints 
 "y >= 0 \n" # add compatibillity constraints 
 "}\n""}\n""}\n""jumps\n""{\n""tran -> rot\n"
 "guard { t = 1  }\n" # edw einai h fasoula gia to jump
 "reset { }\n""parallelotope aggregation { }\n""rot -> tran \n""guard {x = 5}\n""reset{}#x' := x - 4.9 }\n""parallelotope aggregation { }\n""}\n""init\n""{\n""tran \n""{\n"
  "x in [0.0 , 1.0]\n"  # add initial set constraints 
  "y in [0.5 , 1.0]\n"  # add initial set constraints 
  "x3 in [0, 0.1]\n"   # add initial set constraints  
  "}\n""}\n""}\n"
)