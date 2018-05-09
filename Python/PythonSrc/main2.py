#Main2 File: Use Inverting
from DiffusionMaps import *
from GeometricFlows import *
from libfunctions import *
from DataGenerators import *
from animations import *
import sys


Npoints = 400
Nt = 200
tmax = 10
radius = 3
#use a Circle for data 
Data,theta,t = circle_data(Ntheta = Npoints, Nt = 200, tmax = 10)
A = Data.data
GeometricFlow(Data,5000,.005,plot = 1)



# Going to define a heatflow, and save the image
#
#Npoints = 300

#Nt = 200
#Data,theta,t = circle_data(Ntheta = Npoints, Nt = 200, tmax = 10) 


