#Module for different data Declarations
import numpy as np
from DiffusionMaps import mdata
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def circle_data(Ntheta=500,Nt=100,tmax=10,shape = (1,1),plot = 0): #Creates test data for an ellipse using Ntheta uniform samples on the circle and Nt uniform samples in time
    xstretch = shape[0]
    ystretch = shape[1]
    xs = 2*np.pi*np.linspace(1/Ntheta,1,Ntheta) #Even grid from 0 to 2pi
    theta = xs
    #theta = xs - np.sin(xs)/2  # Non-Uniform theta grid
    t = np.linspace(0,tmax,Nt)
    if plot == 1:
        plt.figure()
        plt.plot(xs,theta)
        plt.title(r"Grid Distribution")
        plt.xlabel(r"$\theta$")
        plt.ylabel(r"$f \left (\theta \right )$")
    xn = xstretch*np.cos(theta); yn = ystretch*np.sin(theta);  # Defines data
    xn = xn.reshape(Ntheta,1); yn = yn.reshape(Ntheta,1)
    if plot ==1:
        plt.figure()
        plt.scatter(xn,yn)
        plt.title(r"Data on Manifold")
        plt.xlabel(r"$x_i(\theta)$")
        plt.ylabel(r"$y_i(\theta)$")
    A = np.concatenate((xn,yn),axis = 1)
    X = mdata(A,theta)
    return X,theta,t


def sphere_initial_conditions(Npoints):
    fo = np.ones(Npoints**2)
    fo[100::] = 0
    return fo


def sphere_data(Npoints, Nt, tmax, radius=1,plot = 0,initial= None): #Optional: provide intial condition to plot
    t = np.linspace(0,tmax,Nt)
    xx = np.linspace(0,1,Npoints**2)
    phi = np.pi*np.linspace(1/Npoints,1,Npoints) #Even grid from 0 to 2pi
    theta = 2*np.pi*np.linspace(1/Npoints,1,Npoints) #Even grid from 0 to 2pi
    PHI, THETA = np.meshgrid(phi,theta)
    x = radius*np.sin(PHI)*np.cos(THETA)
    y = radius*np.sin(PHI)*np.sin(THETA)
    z = radius*np.cos(PHI)
    x = x.reshape((Npoints**2,1))
    y = y.reshape((Npoints**2,1))
    z = z.reshape((Npoints**2,1))
    A = np.concatenate((x,y,z),axis=1)
    if plot ==1:
        fig = plt.figure()
        ax = fig.gca(projection = '3d')
        ax.scatter(x,y,z,c = initial,cmap = "coolwarm", vmin = initial.min(), vmax = initial.max())
        plt.title(r"Initial Sphere Data")
        ax.set_xlabel(r'$X$')
        ax.set_ylabel(r'$Y$')
        ax.set_zlabel(r'$Z$')
    plt.show()
    X = mdata(A)
    return X,t      
