# Following code implements the Vicsek model in 2D

# Importing packages
import numpy as np
import matplotlib.pyplot as plt
#import aggregation as agg
import mpl_toolkits.mplot3d.axes3d as axes3d
from scipy.interpolate import griddata
import gudhi
from operator import itemgetter
#import aggregation as agg

class Vicsek2D:
    def __init__(self, N, eta):
        # Initialize simulation
        self.L = 1  # length of the square 2D region to be simulated
        self.halfL = self.L / 2  # half of length (used later for PBCs)
        self.N = N  # number of particles in the 2D region
        self.rho = N/self.L**2  # density of particles in the simulation
        self.eta = e  # noise in the system
        self.r = 0.05  # interaction radius
        self.rsq = self.r * self.r  # square of interaction radius
        self.dt = 1.0  # time step
        self.v = 0.03  # magnitude of velocity
        self.pos = np.random.rand(self.N, 2) * self.L  # random initial position in 2D region
        self.theta = (np.random.rand(self.N) * 2 - 1) * np.pi  # random velocity angle [-pi pi]
        self.vel = np.zeros((self.N, 2))  # initialize velocity array
        self.vel[:, 0] = self.v * np.cos(self.theta)  # velocity along x
        self.vel[:, 1] = self.v * np.sin(self.theta)  # velocity along y
        self.tt = 500 # total number of time steps
        self.rparts = np.eye(N, dtype=np.bool_)  # matrix representing particles within distance r
        
    def main(self):
        plt.ion()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.show(block=False)
        fig.tight_layout()
        fig.canvas.draw()
        plt.show()
        
        
        for t in range(1,self.tt+1,1):
            posn=[]
            
            
            ax.clear()
            #ax.quiver(self.pos[:, 0], self.pos[:, 1], self.vel[:, 0], self.vel[:, 1])
            ax.scatter(self.pos[:, 0], self.pos[:, 1],c='black',s=1)
            #ax.axis(axrange)
            #ax.set_aspect('equal', 'box')
            
            plt.pause(0.01)
            fig.canvas.draw()
            fig.savefig("agg_n_%i_eta_%1.3f_sc_1_t_%i_.png"%(N,e,t) )
            #plt.close()
            #print('vel',self.vel[:,0], self. vel[:,1])
            
            
            for p in self.pos:
                    pt=(p[0], p[1])
                    posn.append(pt)
                #print(posn)    
            post= np.asarray(posn)
            np.savetxt('new_matrix N_%i_t_%i_noise_%1.2f.txt' %(N,t,e),post)     
            self.update()
        

    def update(self):

        # Generate the set of random movements dTheta from [-eta/2, eta/2]
        noise = (np.random.rand(self.N) - 0.5) * self.eta

        # Find particles within distance r
        self.find_particles()

        # Initialize average theta and order
        avg_theta = np.zeros(self.N)
        order=0.0
        for i in range(N):

            # Angles of particles within separation r
            rtheta = self.theta[self.rparts[i, :]]
            avg_theta[i] = np.arctan2(np.mean(np.sin(rtheta)), np.mean(np.cos(rtheta)))

        # Updated angles = avg. angles + noise
        self.theta = avg_theta + noise

        # Updated velocities
        self.vel[:, 0] = self.v * np.cos(self.theta)
        self.vel[:, 1] = self.v * np.sin(self.theta)

        # Updated positions
        self.pos = self.pos + self.vel * self.dt
        
        # order calculation
        

        # Applying periodic boundaries
        #self.pos = np.mod(self.pos, self.L)

    def find_particles(self):

        # Reset rparts matrix
        self.rparts = np.eye(self.N, dtype=np.bool_)

        for i in range(N):
            for j in range(i + 1, N):

                diff = self.pos[i, :] - self.pos[j, :]

                # Apply minimum image criteria for periodic boundaries on diff
                #for dim in range(2):
                    #while diff[dim] > self.halfL:
                        #diff[dim] = diff[dim] - self.L
                    #while diff[dim] < -self.halfL:
                        #diff[dim] = diff[dim] + self.L

                # Separation between particles squared
                sepsq = np.dot(diff, diff)
                rflag = sepsq < self.rsq
                self.rparts[i, j] = rflag
                self.rparts[j, i] = rflag
                

            
        

# Ask user for parameters
N = 500
eta_i= 0.9
eta_f= 1.1


eta = np.arange(eta_i, eta_f,1.0)
for e in eta :
    v2d = Vicsek2D(N,e)

    v2d.main()



