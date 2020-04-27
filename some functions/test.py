import numpy as np 
from numpy import array

from pypoman import compute_polytope_halfspaces
import polytope as pc


vertices = map(array, [1, 3, 2])
A, b = compute_polytope_halfspaces(vertices)

print (A,b)



import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection



P=  [[12.31052917, 45.18707857,  0.549995  ],
        [12.31052917, 45.18707857,  0.550005  ],
        [12.31052917, 45.192353,    0.549995  ],
        [12.31052917, 45.192353 ,   0.550005  ],
        [12.313735  , 45.18707857,  0.549995  ],
        [12.313735  , 45.18707857 , 0.550005  ],
        [12.313735  , 45.192353 ,   0.549995  ],
        [12.313735  , 45.192353 ,   0.550005  ]]  



print (x)
print (y)
print (z)

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from matplotlib import cm


x =[float(item[0]) for item in P]
y =[float(item[1]) for item in P]
z = [float(item[2]) for item in P]
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.2)
plt.show()