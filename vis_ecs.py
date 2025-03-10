#CONSTRUCTING ECS FROM POINT MATRIX
# Importing packages
import numpy as np
import matplotlib.pyplot as plt
#import aggregation as agg
import mpl_toolkits.mplot3d.axes3d as axes3d
from scipy.interpolate import griddata
import gudhi
#import aggregation as agg


def ecs(e):
        # Setup plots
        print('gone, NOISE=',e)
        #plt.ion()
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #plt.show(block=False)
        #fig.tight_layout()
        #fig.canvas.draw()
        #plt.show()
        #axrange = [0, self.L, 0, self.L]
        
        Ms=[]
        Ns = []
        chis =[] 
        bet0=[]
        bet1=[]
        bet2=[]
        sc=6
        tt=500
        N=500
        n1=12
        n2=250
        
        ts=[]
        scale=np.arange(1,sc+1,0.5)
        for t in range(1,tt,2):
            posn=[]
            ts.append(t)
            
            
            for ms in  scale:
                print('scale=',ms, 'time=',t ,'::')
                r_sc= ms/(100)
                b=np.zeros(N-1)
               
                #print(posn)    
                post=np.loadtxt('matrix N_200_t_%i_noise_%1.2f.txt'%(t,e),dtype=float)
                alpha = gudhi.AlphaComplex(post)
                st = alpha.create_simplex_tree(max_alpha_square = r_sc**2)
                simplices = st.get_simplices()

                simp = [s[0] for s in simplices]
                #print (simp)
                #print(sorted(simp))

                #tri = []
                #fig, ax = plt.subplots(1,1)
                #for s in simp:
                    #if len(s)==1: continue
                    #if len(s)==2:
                        #ax.plot(post[s,0], post[s,1], 'r')
                   # else:
                        #tri = plt.Polygon(post[s], fill=True, alpha=0.5)
                       # ax.add_patch(tri)


                #xs, ys = zip(*post)
                #ax.scatter(xs, ys, color='k',s=1)
                #plt.savefig('alpha_vis_eta_%1.2f_t_%i_sc_%1.3f.png' %(e,t,r_sc))
               # plt.close()
                #for i in range(len(points)):
                    #x, y = points[i]
                    #ax.text(x,y,'%i' %i)
                    #plt.show() 
                n=N
                for i in range(0,n-1):     
                    for pt in simp:
                        if len(pt)==i+1:
                            b[i]=b[i]+1

            
                Ms.append(r_sc)
                Ns.append(t)
                BarCodes_1 = st.compute_persistence()
                bet= st.betti_numbers()
                print('betti____',bet)
                b_0= bet[0]
                if (len(bet)>1):
                   b_1= bet[1]
                else:
                   b_1=0
                
               
                bet0.append(b_0)
                bet1.append(b_1)
                #bet2.append(b_2)
                chi_sc=0
                for i in range (0,n-1) :
                    if i%2==0:
                        chi_sc=chi_sc+b[i]
                    else:
                        chi_sc=chi_sc-b[i] 
                
                
                chis.append(chi_sc)
                #bet0.append(b[0])
                #bet1.append(b[1])
                #bet2.append(b[2])
                
                print("**EULER NUMBER**",b[0],b[1],b[2],chi_sc,'ms=',ms,'t=',t)
            
            #N_chi = np.asarray([mS, Ns, chis]).T
        #Beta = np.asarray([ms,Ns, B_0, B_1, chis]).T
        Beta = np.asarray([Ms,Ns, bet0, bet1, chis]).T
        Beta2=np.asarray([Ms,Ns,chis]).T
        
        np.savetxt('component N_%i_sf_%1.3f_noise_%1.2f.txt' %(N,sc,e),Beta)
        np.savetxt('alpha_surface_particle_%i_sf_%1.3f_noise_%1.2f.txt' %(N,sc,e),Beta2)
        
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        xi = np.linspace(min(Ms), max(Ms), n1)
        yi = np.linspace(min(Ns), max(Ns),n2)
        x_arr.append((xi))
        y_arr.append((yi))
        xi,yi = np.meshgrid(xi,yi)

        zi = griddata((Ms,Ns),chis,(xi,yi),method='linear')
        
        

        z_arr.append((zi))
        
       
        ax.plot_surface(xi,yi,zi)

        #ax.scatter(x, y, z,color='red')#,label='data')

        ax.set_xlabel('$r$',size=20)
        ax.set_ylabel('$t$',size=20)
        ax.set_zlabel('$\chi$',size=20)

        plt.savefig('alpha_particle_%i_sf_%1.3f_noise_%1.2f.png' %(N,sc,e))
        #plt.show()
        plt.close()
        

if __name__=='__main__':
   eta_i= 2.1
   eta_f= 5.1
   d_eta=0.1
#grid=900
   x_arr=[]
   y_arr=[]
   z_arr=[]
   eta = np.arange(eta_i, eta_f, d_eta)
   for e in eta :
      ecs(e)

