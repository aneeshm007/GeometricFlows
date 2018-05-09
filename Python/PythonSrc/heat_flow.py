#This is going to construct/plot the heat flow, rather than geometric flow on a manifold
from DiffusionMaps import *
import numpy as np
import matplotlib.pyplot as plt


#A is an mdata class

#let t be a vector 
def heatflow(A,t,f0):

    sol = A.NonUniform_spec_solve_vec(t,f0,eps = .003):

    for i in range(len(t)):
    plt.scatter(A)
