import numpy as np
import gudhi
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
from scipy.interpolate import griddata 

def generate_points(seed, k=4.0, N=10000, ms=1):
        points = [seed]
        #r = ms/(144*4.8)
        for i in range(N):
             x = seed[0] + k*seed[1]*(1-seed[1])
             x = x%1
             y = seed[1] + k*seed[0]*(1-seed[0])
             y = y%1
             p = (x,y)
             #if p not in points:
             points.append(p)
            
             seed = p
    #xs, ys = zip(*points)
    #plt.scatter(xs,ys,s=1)
    #plt.title('k = %2.3f' %k)
    #plt.savefig('acb_%2.3f.png' %k)
    #plt.show()
    #plt.close()
        #points = [p + tuple([r]) for p in points]
        return points
        
ki,kf,dk=4.1,4.2,0.10
si,sf,dsc=1.0/(144*4.8),11/(144*4.8),1.0/(144*4.8)
ks = np.arange(ki, kf, dk)
print(ki,kf,dk)

seed=(0.9, 0.2)
N=10000
dN=500 
ms=1 
#grid=900
x_arr=[]
y_arr=[]
z_arr=[]
ii=0
ii_max=int((kf-ki)/dk)-1
for k in ks:
            
        pointall = generate_points(seed, k, 10000, ms)
        print('k=',k)
        scale=np.arange(si, sf, dsc)
        ms=[]
        Ns = []
        chis =[] 
        B_0=[]
        B_1=[]
        B_2=[] 
        for sc in scale:
            print(sc)       
                                                          
            for n in range(1, N+1, dN):
                print(n)
                b=np.zeros(n-1)
                points = pointall[:n]
                points= np.asarray(points)
                alpha = gudhi.AlphaComplex(points)
                st = alpha.create_simplex_tree(max_alpha_square = sc**2)
                simplices = st.get_simplices()

                simp = [s[0] for s in simplices]
                #print(sorted(simp))

                tri = []
                fig, ax = plt.subplots(1,1)
                for s in simp:
                    if len(s)==1: continue
                    if len(s)==2:
                        ax.plot(points[s,0], points[s,1], 'k')
                    else:
                        tri = plt.Polygon(points[s], fill=True, alpha=0.5)
                        ax.add_patch(tri)


                xs, ys = zip(*points)
                ax.set_xlim(0,1)
                ax.set_ylim(0,1)
                ax.set_box_aspect(1)
                ax.scatter(xs, ys, color='r')
                plt.savefig('alpha_logistics_k_%1.2f_sc_%1.3f_n_%i.png' %(k,sc,n))
                plt.close()
                #for i in range(len(points)):
                    #x, y = points[i]
                    #ax.text(x,y,'%i' %i)
                    #plt.show() 
                for i in range(0,n-1):     
                    for pt in simp:
                        if len(pt)==i+1:
                            b[i]=b[i]+1

            
                ms.append(sc)
                Ns.append(n)
                #BarCodes_Rips1 = st.compute_persistence()
                #bet= st.betti_numbers()
                #print('betti____',bet)
                #b0= bet[0]
                #b1= bet[1]
                #b2=bet[2]
                #B_0.append(b0)
                #B_1.append(b1)
                #B_2.append(b2)
                chi_sc=0
                for i in range (0,n-1) :
                    if i%2==0:
                        chi_sc=chi_sc+b[i]
                    else:
                        chi_sc=chi_sc-b[i] 
                
                
                chis.append(chi_sc)
                
                print("**EULER NUMBER**",chi_sc)    
            #N_chi = np.asarray([mS, Ns, chis]).T
        #Beta = np.asarray([ms,Ns, B_0, B_1, chis]).T
        Beta2=np.asarray([ms,Ns,chis]).T
        #np.savetxt('rips(betti)n_%i_sf_%1.3f_k_%1.2f_ECS.txt' %(N,sf,k),Beta)
        np.savetxt('alpha_surface_n_%i_sf_%1.3f_k_%1.2f.txt' %(N,sf,k) ,Beta2)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        xi = np.linspace(min(ms), max(ms), 300)
        yi = np.linspace(min(Ns), max(Ns), 300)
        x_arr.append((xi))
        y_arr.append((yi))
        xi,yi = np.meshgrid(xi,yi)

        zi = griddata((ms,Ns),chis,(xi,yi),method='linear')
        
        

        z_arr.append((zi))
        
       
        ax.plot_surface(xi,yi,zi)

        #ax.scatter(x, y, z,color='red')#,label='data')

        ax.set_xlabel('$r$',size=20)
        ax.set_ylabel('$t$',size=20)
        ax.set_zlabel('$\chi$',size=20)

        plt.savefig('alpha_n_%i_sf_%1.3f_k_%1.2f.png' %(N,sf,k))
        #plt.show()
        plt.close()
                    #for i in range(len(points)):
                        #x, y = points[i]
                        #ax.text(x,y,'%i' %i)
                        


    
