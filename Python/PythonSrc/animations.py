#Animations
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation, rc

class AnimatedScatter(): #Class to produce an Animated Scatter Plot in 2D or 3D
    """Animates a Scatter Plot/creates an animated heat map given a solution matrix"""
    def __init__(self, sol,A): #
        self.N = len(A) #Number of points
        self.numframes = sol.shape[1]
        self.sol = sol
        self.min = sol.min()
        self.max = sol.max()
        self.xn = A[:,0]
        self.yn = A[:,1]
        dim = A.shape[1]
        self.dim = dim
        if dim == 2:
            self.zn = None
            self.fig, self.ax = plt.subplots()
            plt.xlabel(r"$\theta$")
            plt.ylabel(r"$f\left(\theta \right)$")
            plt.grid()
        elif dim == 3:
            self.zn = A[:,2]
            self.fig = plt.figure()
            self.ax  = self.fig.gca(projection = "3d")
            plt.grid()
        # Then setup FuncAnimation.
        if dim == 3:
            self.ani = animation.FuncAnimation(self.fig, self.update3, init_func = self.setup_plot3,interval=5, blit = False,frames=self.numframes)
        elif dim == 2:
            self.ani = animation.FuncAnimation(self.fig, self.update2, init_func = self.setup_plot2,interval=5, blit = False,frames=self.numframes)

    def setup_plot3(self):
        """Initial drawing of the scatter plot."""

        sol = self.sol
        self.scat = self.ax.scatter(self.xn,self.yn,self.zn, c = sol[:,0], cmap = "coolwarm",vmin = self.min, vmax = self.max)
        self.ax.set_xlabel(r'X')
        self.ax.set_ylabel(r'Y')
        self.ax.set_zlabel(r'Z')
        plt.title("Heatflow in $\mathbb{R}^3$")
        return self.scat,

    def update3(self, i):
        sol = self.sol
        plt.cla()
        self.ax.set_xlabel(r'X')
        self.ax.set_ylabel(r'Y')
        self.ax.set_zlabel(r'Z')
        plt.title("Heatflow in $\mathbb{R}^3$")
        self.scat = self.ax.scatter(self.xn,self.yn,self.zn, c = sol[:,i], cmap = "coolwarm",vmin = self.min, vmax = self.max)
        return self.scat,
    
    def setup_plot2(self):
        sol = self.sol
        self.scat = self.ax.scatter(self.xn,self.yn, c = sol[:,0], cmap = "coolwarm",vmin = self.min, vmax = self.max)

    def update2(self, i):
        sol = self.sol
        plt.cla()
        self.scat = self.ax.scatter(self.xn,self.yn, c = sol[:,i], cmap = "coolwarm",vmin = self.min, vmax = self.max)
        return self.scat,

    def save(self,name,fps,dpi):
        anim2 = self.ani
        anim2.save(name, fps=fps, extra_args=['-vcodec', 'libx264'],dpi = dpi)

    def show(self):
        plt.show()
