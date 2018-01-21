#Library functions
import numpy as np
from sklearn.neighbors import NearestNeighbors
from DiffusionMaps import *

def bandwidthFinder(Data,p,k = False,minkowski_parameter = 2):
    """A function that uses k nearest neighbors to find an epsilon to construct a graph Laplacian:
   Data is either a matrix of points (m x n) ([x_00,x_10,..,x_m0].....[x_0n,x_1n,....x_mn] or
   k is the number of nearest neighbors, default set to ceil(log(N))
   p is the number of points to sample randomly around the manifold
   The optional minkowski parameter is for other metrics of distance
   """
    A = Data.data

    N = len(A)

    if  k:
        k = k

    else:
        k = int(np.ceil(np.log(N)))

    ind = np.random.randint(N,size = p) #This will choose 'p' random indices
    choosePoints = A[ind] #These are the points to choose from.
    neigh = NearestNeighbors(n_neighbors = k, metric ='minkowski', p = minkowski_parameter).fit(A)
    dist, ind = neigh.kneighbors(choosePoints) #Returns the distances
    distance_sum = np.sum(dist)/k #returns the average distances
    epsilon = np.sum(distance_sum)/p
    return epsilon

def polygon_normalization(data): #Can provide a matrix or an mdata object consisting of datapoints in R^n
    """Uses the shoelace formula to find an approximation of the area of the 2-D manifold"""
    if type(data) == mdata:
        A = data.data
    elif type(data) == np.ndarray:
        A = data
    else:
        input_type = str(type(data))
        raise TypeError("Incorrect data type: %s"%input_type)
    #Area Approximation
    Atemp = np.empty((A.shape[0] + 1, A.shape[1])) #Setup for shoelace formula
    Atemp[0:-1,:] = A
    Atemp[-1,:] = A[0]
    A = Atemp
    xadd = A[0:-1,0]; yadd = A[1::,1];
    xsub = A[1::,0]; ysub = A[0:-1,1];
    area = xadd.dot(yadd) - xsub.dot(ysub) # dot products
    return area/2
