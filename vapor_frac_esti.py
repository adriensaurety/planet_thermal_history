import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

G= 6.67430E-11 #gravitational constant
density_planet=5513#kg.m-3


#### fit des hugoniots, here we stock the point and fit we need for the rest of the analysis 

entropy_chondrule=[2028.94,2929.99,3545.65]
temp_chondrule=[2000,4000,6000]
aec,bec=np.polyfit([2028.94,2929.99,3545.65],np.log([2000,4000,6000]),1) #fit give log(temperature) as function of entropy

aep,bep=np.polyfit([4204.95,2367.93],np.log([6000,2000]),1)
apy,bpy=np.polyfit([2000,3000,4000,5000,8000],[63.4,78.3,97.3,120.6,180],1)

P_chond_r0_5_T0_10K=[12.96,
85.8,
108.55,
127.5,
153.14,
182.09,
234.68]
u_chond_r0_5_T0_10K=[0,
2808.44109,
3211.369016,
3550.939248,
4097.908681,
4463.471366,
5185.334982]

a_chond_Pu,b_chond_Pu,c_chond_Pu=np.polyfit(u_chond_r0_5_T0_10K,P_chond_r0_5_T0_10K,2) #this fit give the pressure as a function of the particle velocity

ach,bch=np.polyfit([2000,3000,4000,5000,6000,8000],[85.8,108.55,127.5,153.14,182.09,234.68],1) #this fit give temperature as a function of pressure

print('chondrule hg PT  a=',ach,'b=',bch)
print('chondrule hg ST a=',aec,'b=',bec)
print('chondrule hg Pu a=',a_chond_Pu,'b=',b_chond_Pu,'c=',c_chond_Pu)
#### function

def vesc2v(M,m,v):
    """
    return the v
    elocity of the impactor in m.s-1
    """
    r1=((M/((4/3)*3.1416*density_planet))**(1/3))
    r2=((m/((4/3)*3.1416*density_planet))**(1/3))
    vesc=((2*G*(M+m))/(r1+r2))**(1/2)
    return v*vesc

def u_lim_vap(entropy_vap_lim,aST,bST,aPT,bPT,aPu,bPu,cPu):
    T_lim=np.exp(entropy_vap_lim*aST+bST)
    print(T_lim)
    P_lim=T_lim*aPT+bPT
    print(P_lim)
    u_imp=0
    P_imp=0
    plt.figure()
    while P_imp<P_lim:
        u_imp+=100
        X=np.linspace(0,20000,10000)
        Y=aPu*X**2+bPu*X+cPu
        X_i=u_imp-X
        plt.plot(X,Y)
        plt.plot(X_i,Y)
        for k in range(len(X)):
            if X[k]<X_i[k]+3.0 and X[k]>X_i[k]-3.0:
                P_imp=Y[k]
    return u_imp
#### read the file

data=np.loadtxt('sim_col_nojupiter.dat', usecols=[0,3,4,12],skiprows=1)

data[:,0]=data[:,0]*24*3600*6/365.25
data[:,1]=data[:,1]*2E30
data[:,2]=data[:,2]*2E30

for k in range(len(data[:,3])):
    data[k,3]=vesc2v(data[k,1],data[k,2],data[k,3])


#### criteria to assume vaporization
entropy_vap_lim=3200

u_lim_chond=u_lim_vap(entropy_vap_lim,aec,bec,ach,bch,a_chond_Pu,b_chond_Pu,c_chond_Pu)
print(u_lim_chond)
u_lim_pyr=u_lim_vap(entropy_vap_lim,aep,bep,apy,bpy,a_chond_Pu,b_chond_Pu,c_chond_Pu)
print(u_lim_pyr)

#### loop over all the impacts
ilv_chond=0
for imp in data[:,3]:
    if imp>u_lim_chond:
        ilv_chond+=1
ilv_pyr=0
for imp in data[:,3]:
    if imp>u_lim_pyr:
        ilv_pyr+=1

print("Chondrule case: proportion of impact leading to vaporization=",ilv_chond/len(data[:,3]))

print("Pyrolite case: proportion of impact leading to vaporization=",ilv_pyr/len(data[:,3]))
#### plot

plt.figure()
N,bins,patches=plt.hist(data[:,3], 100,ec="k")
cmap=plt.get_cmap('jet')
none="#6A5D7B"
chondrule="#749C75"
chond_pyr='#E9D985'
for i in range(len(bins)-1):
    if bins[i]<u_lim_pyr:
        patches[i].set_facecolor(none)
    if bins[i]>u_lim_pyr and bins[i]<u_lim_chond:
        patches[i].set_facecolor(chond_pyr)
    if bins[i]>u_lim_chond:
        patches[i].set_facecolor(chondrule)

handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in [none,chond_pyr,chondrule]]
labels= ["no vaporization","partial vap of pyrolite","partial vap of chondrule and pyrolite"]
plt.legend(handles, labels)

plt.xlabel("impactor velocity (m/s)", fontsize=16)  
plt.ylabel("number of impact", fontsize=16)
plt.text(15000,25,"impact% leading to vaporization if \nall bodies are chondrules="+str(100*round(ilv_chond/len(data[:,3]),4))+"%")
plt.xticks(fontsize=14)  
plt.text(15000,18,"impact% leading to vaporization if \nall bodies are pyrolite="+str(100*round(ilv_pyr/len(data[:,3]),3))+"%")
plt.yticks(fontsize=14)
plt.savefig('histo_velocity_vap.png')
plt.show()
