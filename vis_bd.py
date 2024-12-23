# Following code implements the Vicsek model in 2D from PRL 75 1226 (1995)
# For possible use in AM 115

# Importing packages
import numpy as np
import matplotlib.pyplot as plt
#import aggregation as agg
import mpl_toolkits.mplot3d.axes3d as axes3d
from scipy.interpolate import griddata
import gudhi
from gudhi.wasserstein import wasserstein_distance
from operator import itemgetter
import math
from numpy import inf
#import aggregation as agg

def wd(e,t):
        # Setup plots
        print('gone','e=',e)
        #plt.ion()
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #plt.show(block=False)
        #fig.tight_layout()
        #fig.canvas.draw()
        #plt.show()
        #axrange = [0, self.L, 0, self.L]
        
      
        sc=80
        #for ms in range(self.sc,self.sc+1,1):
        for tt in range(t,t+1,1):
        
            for ms in range(sc,sc+1,1):
            
                post=np.loadtxt('matrix N_500_t_%i_noise_%1.2f.txt'%(tt,e),dtype=float)
                #ax.clear()
                #ax.quiver(self.pos[:, 0], self.pos[:, 1], self.vel[:, 0], self.vel[:, 1])
                #ax.scatter(self.pos[:, 0], self.pos[:, 1],c='black',s=ms)
                #ax.axis(axrange)
                #ax.set_aspect('equal', 'box')
            
                #plt.pause(0.01)
                #fig.canvas.draw()
                #fig.savefig("agg_n_%i_eta_%1.3f_sc_%i_t_%i_.png"%(N,e,ms,t) )
                #plt.close()
                r_sc = ms/(100) ###
                
                #for p in self.pos:
                    #pt=(p[0], p[1])
                    #posn.append(pt)
                #print(posn)    
                #post= np.asarray(posn)
                alpha = gudhi.AlphaComplex(post)
                st = alpha.create_simplex_tree(max_alpha_square = r_sc**2)
                simplices = st.get_simplices()

                simp = [s[0] for s in simplices]
                #print (simp)
                #print(sorted(simp))

                #tri = []
                fig, ax = plt.subplots(1,1)
                #for s in simp:
                   # if len(s)==1: continue
                    #if len(s)==2:
                        #ax.plot(post[s,0], post[s,1], 'r')
                    #else:
                        #tri = plt.Polygon(post[s], fill=True, alpha=0.5)
                        #ax.add_patch(tri)


                xs, ys = zip(*post)
                ax.scatter(xs, ys, color='k',s=ms)
                #plt.show()
                #plt.savefig('alpha_vis_eta_%1.2f_sc_%1.3f_n_%i.png' %(e,r,t))
                plt.close()
                pdn=st.persistence()
                #pdn2=st2.persistence()
                #print(pdn)  
            
                if(len(pdn)>0):
                    print('pdn loop')
                    #print(pdn)
                    gudhi.plot_persistence_diagram(np.srt.(pdn))
                    plt.savefig('pd_N=500_t=%i_eta_%1.2f.png'%(t,e))
                    gudhi.plot_persistence_barcode(np.sqrt(pdn))
                    plt.savefig('bc_N=500_t=%i_eta_%1.2f.png'%(t,e))
                    #plt.close()
                    pd_intervals1 = np.sqrt(st.persistence_intervals_in_dimension(1))
                    pd_intervals0 = np.sqrt(st.persistence_intervals_in_dimension(0))
                    pdl1.append((pd_intervals1))
                    pdl0.append((pd_intervals0))
                    


if __name__=='__main__':
   eta_i=0.0
   eta_f=0.95
   d_eta=0.9
   tt=np.arange(1,500,2)
   #tf=301
   bdt=0.0
   wt=0.0
   eta = np.arange(eta_i, eta_f, d_eta)
#tt= np.arange(ti,tf,dt)
   #pdl0=[]
   #pdl1=[]
   for t in tt :
    #wd1=[]
    #wd0=[]
         pdl0=[]
         pdl1=[]
         for e in eta:
        
            wd(e,t)
    
#print(pdl1, pdl0)
         print("t==",t, len(pdl1), len (pdl0))            
         pd11_intervals = np.asarray(pdl1[0])

         pd12_intervals = np.asarray(pdl1[1])
    
         wasserstein_distance1 = wasserstein_distance(pd11_intervals, pd12_intervals,order=1,internal_p=2.)
         bd1= gudhi.bottleneck_distance(pd11_intervals, pd12_intervals)  
         pd01_intervals = np.asarray(pdl0[0])

         pd02_intervals = np.asarray(pdl0[1])
#print ('before pd01',pd01_intervals)
#print ('before pd02',pd02_intervals)
#pd01_intervals[pd01_intervals[:,1]==inf,1]=1  
#pd02_intervals[pd02_intervals[:,1]==inf,1]=1 
#print('after pd01',pd01_intervals) 
#print('after pd02',pd02_intervals) 
         max01= pd01_intervals.max()
         max02= pd02_intervals.max()
         max11= pd11_intervals.max()
         max12= pd12_intervals.max()
         min01= pd01_intervals.min()
         min02= pd02_intervals.min()
         min11= pd11_intervals.min()
         min12= pd12_intervals.min()
         maxf= max([ max11, max12] )
         minf= min([ min11, min12] )
         print('max==', max01, max02,max11, max12, maxf) 
         print('min==', min01, min02,min11, min12, minf) 
         wasserstein_distance0=wasserstein_distance(pd01_intervals, pd02_intervals,order=1,internal_p=2.)
         bd0=gudhi.bottleneck_distance(pd01_intervals, pd02_intervals)
         WD=np.asarray([e,t,wasserstein_distance0,wasserstein_distance1,bd0,bd1]).T
         print('wasserstein distance=',wasserstein_distance0,wasserstein_distance1,'t=',t)
         np.savetxt('WD_time_%i_noise_%1.2f_rad_0.05.txt' %(t,e),WD) 
         len_bd0= max(len(pd01_intervals),len(pd02_intervals))
         len_bd1= max(len(pd11_intervals),len(pd12_intervals))
         sumw=(2*(wasserstein_distance0+wasserstein_distance1 ))
         print('sum_wass',sumw,'t=',t)  
         sumbd=math.sqrt(2)*( len_bd0* math.sqrt(bd0) + len_bd1* math.sqrt(bd1))
    #plt.close()

         print('sum_bd',sumbd,'t=',t)
         with open('ph_t_%i'%(t),'w') as f:
            f.write(str(sumw))
            f.write('\t')
            f.write(str(sumbd))
            f.write('\t')
            f.write(str(minf))
            f.write('\t')
            f.write(str(maxf))
            
            

         f.close()
         bdt=bdt+sumbd
         wt= wt+sumw
   print('total sum wass, bd', wt, bdt)
           
