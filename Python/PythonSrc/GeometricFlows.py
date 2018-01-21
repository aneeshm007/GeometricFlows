#GeometricFlows
from libfunctions import *
import numpy as np

def GeometricFlow(manifold_data,T,dt,plot = 0): # Runs a Geometric flow for "T" iterations of time step dt
    data = manifold_data
    N = data.length
    solx = np.empty((N,T))
    soly = np.empty((N,T))
    A = data.data
    tempx = A[:,0]; tempy = A[:,1]
    npoints = int(np.ceil(np.sqrt(data.length)))
    initial_area = polygon_normalization(data)
    print(initial_area)
    for i in range(T):
        band_eps = bandwidthFinder(Data = data, p = npoints)
        #if band_eps < 10**-5:
        #    band_eps = 10**-5
        #tempx = data.NonUniform_spec_solve((t[i+1]-t[i]),f0 = tempx,eps = band_eps)
        #tempy = data.NonUniform_spec_solve((t[i+1]-t[i]),f0 = tempy,eps = band_eps)
        tempx = data.NonUniform_spec_solve(dt,f0 = tempx,eps = band_eps)
        tempy = data.NonUniform_spec_solve(dt,f0 = tempy,eps = band_eps)
        area_norm = initial_area/polygon_normalization(data)
        print("norm: ", area_norm)
        print("iteration: ", i)
        print("current_area: ", polygon_normalization(data))
        print("current epsilon: ", band_eps)
        tempx = tempx*(np.sqrt(area_norm))
        tempy = tempy*(np.sqrt(area_norm))
        if plot == 1:
            plt.plot(tempx,tempy)
            plt.savefig('p%d'%i)
        solx[:,i] = tempx
        soly[:,i] = tempy
        temp = np.concatenate((tempx.reshape((N,1)),tempy.reshape((N,1))),axis = 1)
        data.data = temp


    return solx,soly
