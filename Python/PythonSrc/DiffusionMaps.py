#Geometric Flows Python Code
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
import scipy as sp
from matplotlib import animation, rc
from copy import deepcopy

class mdata: # Defines class for a set of data values
    def __init__(self,Data,theta = False, phi = False): #data is numpy array of points that lie on some manifold
        self.data = Data #The data contains some data "sampled from a manifold". Simply consists of points in R^n
        self.theta = theta #This parameter is necessary to plot and only applies to a circle
        self.phi = phi #Phi for plotting the eigenfunctions
        self.identity = np.identity(self.data.shape[0]) #Identity matrix of the size of data
        self.length = len(self.data) #This the number of samples of the manifold we have
        return

    def graph(self):
        return cdist(self.data,self.data) #cdist operates the same as pdist2 in matlab. It returns a matrix whose
                                          #(i,j) entry is the Euclidean distance between points i and j

    def kernel(self,eps):     
        graph = self.graph() # The kernel magnifies the distances in the weighted graph
        return np.exp(np.square(graph)/(-2)/eps/eps)

    def dmatrix(self,eps): #The diagonal entries in the dmatrix correspond to the row sums of the kernel
        return np.diag(np.sum(self.kernel(eps),axis=1))

    def laplacian(self,eps): #Returns the unnormalized graph Laplacian
        D = self.dmatrix(eps)
        K = self.kernel(eps)
        L = (D-K)/eps**2
        return L

    def leftright(self,eps): #Does a normalization for non-uniform sampling
        x = np.shape(self.data)
        Ldim = (x[0],x[0])
        I = np.identity(Ldim[0])
        K = self.kernel(eps)
        D = self.dmatrix(eps)
        Dinv = np.linalg.inv(D)
        Ktilde = K.dot(Dinv)
        Dtilde = np.diag(np.sum(Ktilde,axis=1))
        Dtildeinv = np.linalg.inv(Dtilde)
        #Left normalization
        Khat = Dtildeinv.dot(Ktilde)
        Lhat = (I - Khat)/(eps**2)
        return Lhat

    def symmetric(self,eps):
        x = np.shape(self.data) #Returns a symmetric matrix so we can use a faster eigensolver
        Ldim = (x[0],x[0])
        I = np.identity(x[0])
        K = self.kernel(eps)
        D = self.dmatrix(eps)
        Dinv = np.linalg.inv(D)
        Ktilde = K.dot(Dinv)
        Ktilde = Dinv.dot(Ktilde)
        Dtilde = np.diag(np.sum(Ktilde,axis=1))
        Dtildeinv = np.linalg.inv(Dtilde)

        Dhalf = np.power(np.diag(Dtildeinv),0.5)
        Dhalf = np.diag(Dhalf)
        #Symmetric normalization

        Ktilde = Dhalf.dot(Ktilde)
        Ktilde = Ktilde.dot(Dhalf)
        Lsym = I - Ktilde

        Lsym = (Lsym + Lsym.T)/2
        Lsym = Lsym/(eps**2)

        return(Lsym,Dhalf)


    def eigs(self,eps = .003): #Returns the eigenvalues using the symmetric normalization
        emach = np.finfo(float).eps

        Lsym,Dhalf = self.symmetric(eps)
        eigenValues, eigenVectors = np.linalg.eigh(Lsym)
        idx = eigenValues.argsort()
        eigenValues = eigenValues[idx]
        eigenVectors = eigenVectors[:,idx]
        return eigenValues, Dhalf.dot(eigenVectors)

    def eigs2(self,eps=.003): #nonsymmetric eigenvalue
        L = self.leftright(eps)
        eigenValues, eigenVectors = sp.linalg.eig(L)
        idx = eigenValues.argsort()
        eigenValues = eigenValues[idx]
        eigenVectors = eigenVectors[:,idx]
        return np.real(eigenValues), np.real(eigenVectors)

    def eigplot2(self,n,eps): #Plots the eigenfunctions of the Laplacian of a radially parameterized manifold
        if not self.theta:
            return
        eigenValues,eigenVectors = self.eigs(eps)
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.8, hspace=0.3)
        N = eigenVectors.shape[1]
        for i in range(n):
            #plt.subplot(np.ceil((N+1)/2),2,i+1)
            plt.figure()
            plt.grid(1)
            plt.plot(eigenVectors[:,i])
            plt.title(r"%d eigenfunction of L" %i)
            plt.xlabel(r"t")
            plt.ylabel(r"$\Delta_g \vec{f}$")
        return

    def eigplot3(self,n,eps): #Plots the eigenfunctions of the Laplacian of a spherically paramaterized manifold
        if (not self.theta) or (not self.phi):
            return
        eigenValues,eigenVectors = self.eigs(eps)
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.8, hspace=0.3)
        N = eigenVectors.shape[1]
        for i in range(n):
            #plt.subplot(np.ceil((N+1)/2),2,i+1)
            plt.figure()
            plt.grid(1)
            plt.plot(eigenVectors[:,i])
            plt.title(r"%d eigenfunction of L" %i)
            plt.xlabel(r"t")
            plt.ylabel(r"$\Delta_g \vec{f}$")
        return
    def NonUniform_spec_solve_vec(self,t,f0,eps = .003): #Will solve the heat equation for all time vectors in t, and return a matrix
                                                         # "sol" whose ith column is the solution for t[i]. For solving heat equation
        N = self.length
        cols = len(t)
        sol = np.empty([N,cols])
        V,U = self.eigs(eps)
        D = self.dmatrix(eps)
        K = self.kernel(eps)
        Dinv = sp.linalg.inv(D)
        Ktilde = K.dot(Dinv)
        Ktilde = Dinv.dot(Ktilde)
        Dtilde = np.diag(np.sum(Ktilde,axis=1))
        E = ((1/N)*(U.T).dot(Dtilde)).dot(U)
        Einv = sp.linalg.inv(E)
        Einv = np.diag(np.sqrt(np.diag(Einv)))  #The matrix multiplication is done in several steps but it will save some
        Utilde = U.dot(Einv)  #recomputation
        for i in range(cols):
            v1 = np.exp(-1*(t[i])*V)
            v1 = np.diag(v1)
            fhat = (Utilde.T.dot(Dtilde)).dot(f0)/N
            temp1 = Utilde.dot(v1)
            sol[:,i] = temp1.dot(fhat)
        return sol

    def NonUniform_spec_solve(self,t,f0,eps = .003): #Will solve the heat equation for all time vectors in t, and return a matrix
                                                     # "sol" whose ith column is the solution for t[i]
        N = self.length
        V,U = self.eigs(eps)
        D = self.dmatrix(eps)
        K = self.kernel(eps)
        Dinv = sp.linalg.inv(D)
        Ktilde = K.dot(Dinv)
        Ktilde = Dinv.dot(Ktilde)
        Dtilde = np.diag(np.sum(Ktilde,axis=1))
        E = ((1/N)*(U.T).dot(Dtilde)).dot(U)
        Einv = sp.linalg.inv(E)
        Einv = np.diag(np.sqrt(np.diag(Einv)))  #The matrix multiplication is done in several steps but it will save some
        Utilde = U.dot(Einv)  #recomputation
        v1 = np.exp(-1*(t)*V)
        #v1 = np.exp(-1*t*V - t*np.square(V)/(1000**2))
        v1 = np.diag(v1)
        fhat = (Utilde.T.dot(Dtilde)).dot(f0)/N
        temp1 = Utilde.dot(v1)
        sol = temp1.dot(fhat)
        print("max eigenvalue: ", np.max(V))
        return sol
