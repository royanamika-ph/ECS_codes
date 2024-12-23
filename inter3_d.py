#computes EM between two ECSs with k1 and k2
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import math

x1 = []
y1= []
z1 = []
x2 = []
y2= []
z2 = []
x_arr=[]
y_arr=[]
z_arr=[]
n=3000
k1=4.60
k2=4.50
 #first data file 
fr1=open("chi_with_time_k_%1.2f_ECS.txt"%(k1),"r")

while True:
    line=fr1.readline()
    if not line:
        break
    x1.append(float(line.split()[0]))
    y1.append(float(line.split()[1]))
    z1.append(float(line.split()[2]))
fr1.close()

#second data file
fr2=open("chi_with_time_k_%1.2f_ECS.txt"%(k2),"r")
while True:
    line=fr2.readline()
    if not line:
        break
    x2.append(float(line.split()[0]))
    y2.append(float(line.split()[1]))
    z2.append(float(line.split()[2]))
fr2.close()

# surface plot interpolation

xi = np.linspace(min(x1), max(x1), n)
yi = np.linspace(min(y1), max(y1), n)
x_arr.append((xi))
y_arr.append((yi))
xi,yi = np.meshgrid(xi,yi)
zi = griddata((x1,y1),z1,(xi,yi),method='linear')
z_arr.append((zi))

xi2 = np.linspace(min(x2), max(x2), n)
yi2 = np.linspace(min(y2), max(y2), n)
x_arr.append((xi2))
y_arr.append((yi2))
xi2,yi2 = np.meshgrid(xi2,yi2)
zi2=griddata((x2,y2),z2,(xi2,yi2),method='linear')
z_arr.append((zi2))
#xf = np.linspace(min(x2), max(x2), 10) 
#yf = np.linspace(min(y2), max(y2), 10) 
#zf = np.linspace(max(z1),min(z2), 10)

#X, Y, Z = np.meshgrid(x1, y1, z1) 
# plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot_surface(xi,yi,zi,alpha=0.7)
ax.plot_surface(xi2,yi2,zi2,alpha=1)


print('arr',z_arr)


ax.set_xlabel('$r $',fontsize=20, labelpad=7)
ax.set_ylabel('$t $',fontsize=20, labelpad=9)
ax.set_zlabel('$\chi$ ',fontsize=20, labelpad=7)
plt.savefig('replot_double_surface %1.2f_%1.2f.png'%(k1,k2),dpi=100)
#plt.close(fig)
plt.show()

#Using Monte Carlo integration approach
for ical in range (0,1):
        dA=((np.max(x_arr)- np.min(x_arr))*(np.max(y_arr)-np.min(y_arr)))/(n*n)
        print(np.max(x_arr),np.min(x_arr),np.max(y_arr),np.min(y_arr),dA)
        dV1=0.0
        
        count=0
        print('ical=',ical,"len_z=",len(z_arr[ical]),(z_arr[ical].ndim),(z_arr[ical].shape))
        for xx in range(0,len(z_arr[ical])):
            for yy in range(0,len(z_arr[ical])):
                d_chi= (z_arr[ical+1])[xx][yy]-(z_arr[ical])[xx][yy]
                count=count+1
                
                dV1= dV1+ (d_chi**2)*dA #sqr grid vol sum
                #print((z_arr[ical+1])[xx][yy],(z_arr[ical])[xx][yy],'d_chi',d_chi,'dv',(d_chi*dA)**2)
        dV1=math.sqrt(dV1) #L_2 Euler Metric Distance
        
        print("int",dV1,'count=',count)
        with open('newvol_%1.2f-%1.2f.txt'%(k2,k1),'w') as f:
            f.write(str(k2)+'-'+str(k1))
            f.write('\t')
            f.write(str(dV1))

        #np.savetxt('ecs %i.txt' %ical, dV1,dV2) 
        f.close()    
