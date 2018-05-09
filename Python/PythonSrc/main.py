#Main File
from DiffusionMaps import *
from GeometricFlows import *
from libfunctions import *
from DataGenerators import *
from animations import *
import sys

#Just uses a sphere but really you can use any manifold embedded in R3 to get an animation


args = sys.argv
Npoints, Nt, tmax, radius  = args[1::]


print(args)


fo = sphere_initial_conditions(int(Npoints))
Data,t = sphere_data(Npoints=int(Npoints),Nt=int(Nt),tmax=int(tmax),radius =int(radius), plot=0,initial=fo)
A = Data.data

eps_value = bandwidthFinder(Data,p=20,k = False,minkowski_parameter = 2)
#sol = Data.NonUniform_spec_solve_vec(eps = eps_value,t = t,f0 = fo)

#X = AnimatedScatter(sol,A)
#
#X.save("sphere4.mp4",20,200)

#Geometric Flow Example
T = 10
dt = .1    
#solx,soly = GeometricFlow(manifold_data,T,dt)

GeometricFlow(Data,T,dt,plot = 1)

