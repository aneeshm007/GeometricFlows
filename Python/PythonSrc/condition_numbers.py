from BandwidthFinder import *
DiffusionMaps import * from GeometricFlows import * 
from libfunctions import * 
from DataGenerators import * 
from animations import * 
import sys 
import random 

"""
  Idea is compute dtilde(zi,xj), where $xj$ is the exact data set
  
  Then will calculate Ktilde = e^(-dtilde/epsilon), where epsilon is fixed
  Then will calculate cond(Ktilde). Goal: Move z_i in such a way to minimize condition number.
  Weyl's theorem: Small perterbation of matrix results in small pertebations of the eigenvalues
  Gradient Descent: c = CondNum(Z)
  gradc =  (CondNum(Z + epsilonEij) - CondNum(Z - epsilonEij))/(2*epsilon) 
  c_{k+1} = c_k + n*gradc
  
"""





