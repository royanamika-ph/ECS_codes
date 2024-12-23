#Codes used for building ECS with coarse-graining approach
import numpy as np
import matplotlib.pyplot as plt
import aggregation as agg   #the fortran subroutine used for cluster counting and Euler characteristic
import mpl_toolkits.mplot3d.axes3d as axes3d
from scipy.interpolate import griddata
#from matplotlib.tri import Triangulation

def generate_points(seed, k=4.0, N=1000, ms=1):
    points = [seed]
    r = ms/(144*4.8) #converting pixel size to euclidean distance
    for i in range(N):
        x = seed[0] + k*seed[1]*(1-seed[1])
        x = x%1
        y = seed[1] + k*seed[0]*(1-seed[0])
        y = y%1
        p = (x,y)
        points.append(p)
        seed = p
    #xs, ys = zip(*points)
    #plt.scatter(xs,ys,s=1)
    #plt.title('k = %2.3f' %k)
    #plt.savefig('acb_%2.3f.png' %k)
    #plt.show()
    #plt.close()
    points = [p + tuple([r]) for p in points]
    return points

def chi_t_k(seed=(0.9, 0.2), ki=4.1, kf=4.2, dk=0.1, N=10000, dN=500, ms=11, grid=900):
    ''' calculate chi for the system evolving with time and size with a fixed k '''
    ks = np.arange(ki, kf, dk)
    #points = generate_points(seed, k, N, ms)
    x_arr=[]
    y_arr=[]
    z_arr=[]
    ii=0
    ii_max=int((kf-ki)/dk)-1
    for k in ks:
        Ns = []
        chis = []
        mS=[]
        print('k=',k)
        for isize in range (1,ms): 
            print('size=',isize)
            points = generate_points(seed, k, N, isize)
            for n in range(1, N+1, dN):
                pointn = points[:n]
                mat = agg.pixelize(pointn, grid)  #using fortran subroutine to discretize in grids
                plt.imshow(mat,origin='lower', cmap='gray')
                plt.savefig('aggregation_N_%i_ms_%i_k_%1.2f.png' %(n,isize,k))
              
                plt.close()

            
                chi_sq = agg.sq_chi(mat)
                mS.append(isize)
                Ns.append(n)
                chis.append(chi_sq)
        N_chi = np.asarray([mS, Ns, chis]).T
        np.savetxt('chi_with_time_k_%1.2f_ECS.txt' %k ,N_chi)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        xi = np.linspace(min(mS), max(mS), 3000)
        yi = np.linspace(min(Ns), max(Ns), 3000)
        x_arr.append((xi))
        y_arr.append((yi))
        xi,yi = np.meshgrid(xi,yi)

        zi = griddata((mS,Ns),chis,(xi,yi),method='linear')

        print('xi=',xi,'yi=',yi,'zi=',zi)

        z_arr.append((zi))

        ax.plot_surface(xi,yi,zi)

        ax.set_xlabel('$scale \ (r) $',fontsize=12, labelpad=7)
        ax.set_ylabel('$time\ step \ (t) $',fontsize=12, labelpad=9)
        ax.set_zlabel('$Euler \ Characteristics \ (\chi)$ ',fontsize=12, labelpad=7)

        plt.savefig('replotted_chi_vs_N_k_%1.2f.png' %k) #the ECSplotting
       # plt.show()
        plt.close()
    
               
        
    
if __name__ == '__main__':
    ki,kf,dk=4.1,4.2,0.1
    ks = np.arange(ki, kf, dk)
    print(ki,kf,dk)
    chi_t_k(ki=4.1, kf=4.2, dk=0.1, grid=900)
    
    
