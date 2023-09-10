import numpy as np
import matplotlib.pyplot as plt
import math
"""
This code analyzed data of impact in planetary disk using the hugoniots A. Saurety and R Caracas computed.

It was made for data of S. Raymond (Université de Bordeaux). N-body simulation of accretion of proto-planetary disk. (Raymond et al 2009)

Assumption: the bodies are differenciated into core (pure Fe) and mantle (chondrule compostion),  Al26 heating do not exist anymore, the cooldown of a planet is pure radiative cooldown and the temperature of the whole surface after an impact is the peak temperature of the hugoniot.

author: A. Saurety
date: september 
"""

#### constant

kb=1.38E-23
silicate_wt=0.6
density_planet=5513#kg.m-3
molar_mass_mantle=5.1574 #of chondrule kg.mol-1
molar_mass_core=58 #pure Fe
avogadro=6E23
G= 6.67430E-11 #gravitational constant
Cp= 1200 #Heat capacity J.kg-1.K-1 for pyrolite according Sixtrude 2021

#### fit hugoniots, here we stock hugoniot we need for the rest of the analysis P=pressure, T=temperature, u=particle velocity, S=entropy
# those fits have been computed from the data of A. Saurety and R. Caracas

#fit in the Pu digram with P=a*u**2+b*u+c (forbes 2012) (P in GPa, u in m.s-1)

caPu=6.8045e-6
cbPu=0.0073486
ccPu=12.9917   #for chondrule composition (mantle)

#faPu=
#fbPu=
#fcPu=     #for Fe (core)
#
#fit in the PT diagram with T=a*P+b (P in GPa, T in K)

caPT=0.024916
cbPT=32.35  #chondrule

#faPT=
#fbPT=   #Fe

#fit in the ST diagram with ln(T)=a*S+b (T in K, S in J/kg/K)

caST=0.0007276
cbST=6.1354 #chondrule

#faST=
#fbST=    #Fe

#### global function

def extract_SR_data(file):
    """
    extraction of data
    """
    with open(file) as f:
        lines=f.read().splitlines()
        list_name=[]
        for k in range(1,len(lines)): 
            list_name.append(lines[k].split()[1])   #name 

    data=np.loadtxt(file, usecols=[0,3,4,12],skiprows=1)
    data[:,0]=data[:,0]*24*3600*6
    data[:,1]=data[:,1]*2E30
    data[:,2]=data[:,2]*2E30
    return list_name,data


#### class

class planet():
    """
    this class represent the object planet. A planet has a temperature, a mass, a radius, a mass of silicate, a thermal history.

    During an impact, the temperature of the target is updated (it has cooldown since the last impact), the mass is updated (merging of the two body). Then the heat produced is computed. and the new temperature is updated.

    To compute the cooldown we use a radiative cooldown of a sphere (black body model)
    """
    def __init__(self,name,M):
        self.T=0
        self.name=name
        self.impact_time_list=[]
        self.temperature_list=[]
        self.M=M
        self.radius=((M/((4/3)*3.1416*density_planet))**(1/3)) #radius in meter, using the density of the present Earth
        self.silcont=M*silicate_wt
        self.last_imp=0
    def update_mass(self,m): 
        """
        update_mass(m)
        no return
        update the mass of the body, because of accretion it is different at each impact
        need the mass of the impactor m
        """
        M=self.getmass()
        new_mass=M+m
        self.M=new_mass
        self.radius=((new_mass/((4/3)*3.1416*density_planet))**(1/3)) #radius in meter, using the density of the present Earth
        self.silcont=new_mass*silicate_wt
    
    def getname(self):
        return self.name
    
    def getsil(self):
        return self.silcont
    
    def getT(self):
        return self.T
    
    def get_last_imp(self):
        return self.last_imp

    def getradius(self):
        return self.radius

    def getmass(self):
        return self.M
    
    def getimplist(self):
        return self.impact_time_list

    def getTlist(self):
        return self.temperature_list
    
    def vesc2v(self,v_adim,m):
        """
        transform the adimensional impact velocity express relatively to the escape velocity  in a velocity in m.s-1
        """
        M=self.getmass()
        r1=((M/((4/3)*3.1416*density_planet))**(1/3))
        r2=((m/((4/3)*3.1416*density_planet))**(1/3))
        vesc=((2*G*(M+m))/(r1+r2))**(1/2)
        return v_adim*vesc

    def v_to_T(self,v,aPu,bPU,cPU,aPT,bPT):
        P_imp=0
        X=np.linspace(0,80000,20000)
        X_i=v-X
        for k in range(len(X)):
            P_dir=aPu*X[k]**2+bPu*X[k]+cPu
            P_inv=aPu*X_i[k]**2+bPu*X_i[k]+cPu
            if P_inv<P_dir[k]+1.0 and P_inv>P_dir-1.0:
                P_imp=P_dir
        T_imp=P_imp*aPT+bPT
        return T_imp

    def cooldown(self,t):  
        """
        cooldown(t)
        t is the time of the new impact in s 
        this function compute the temperature after a time t
        """
        T=self.getT()
        silcont=self.getsil()
        r=self.getradius()
        M=self.getmass()
        last_imp_time=self.get_last_imp()
        if T!=0:
            self.T=(1/((((t-last_imp_time)*2*1*5.6703E-8*4*3.1416*r**2)/(kb*(silcont*248*avogadro/molar_mass_mantle+(M-silcont)*avogadro/molar_mass_core)))+(1/(T**3))))**(1/3)  #the new temperature depends of the time between two impacts

    def heating_core(self,T_core,v):
        """
        function to compute the heating of the core of the impacted body with a mass M during the shock with an impactor of mass m at velocity v.

        the core as an initial temperature of T_core
        """
    
        return T_core

    def heating_mantle(self,T_mantle,v):
        """
        function to compute the heating of the mantle of the impacted body during the shock with an impactor at velocity v.

        the mantle as an initial temperature of T_mantle
        """
        if T_mantle>1000:
            print('WARNING USE the hot hugoniot')
        
        
         
        return T_mantle


    def impact(self,t,M,m,v_adim):  #M=mass impacted, m=mass impactor
        """
        function that update the mass and the temperature of the planet after an impact occuring at time t between an impacted body with mass M and impactor of mass m at velocity v
        """
        self.cooldown(t)    #computed the Temperature just before the impact
        T=self.getT()        
        self.update_mass(m) #update mass to get the mass at the time of the impact
        implist=self.getimplist()
        Tlist=self.getTlist()
        self.impact_time_list.append(t/(365*24*3600))
        self.temperature_list.append(T)
        M=self.getmass()
        v=self.vesc2v(v_adim,m)
        energy_kin=0.5*m*v**2  #has to be in joules
        new_T=T+energy_kin/(((m+M)*Cp))  #the new temperature is the initial+ the heat du to the impact energy
        self.T=new_T
        self.last_imp=t+1
        self.impact_time_list.append(self.get_last_imp()/(365*24*3600))
        self.temperature_list.append(self.getT())

    def plot_thermal_history(self):
        implist=self.getimplist()
        Tlist=self.getTlist()
        plt.figure()
        plt.plot(implist,Tlist,label=str(self.getname()))
        plt.xlabel('time (years)')
        plt.ylabel('temperature (K)')
        plt.legend()



####
if __name__=='__main__':
    list_name,data=extract_SR_data('4/sim_allcollisions.dat')
    print(list_name)
    print(len(list_name))
    bodies=[]
    nb_impact=len(list_name)
    bodies.append(planet(list_name[0],data[0,1]))
    for k in range(1,nb_impact):
        if list_name[k]!=list_name[k-1]:
            bodies.append(planet(list_name[k],data[k,1]))

    for time in range(len(list_name)):
        c=0
        while bodies[c].getname()!=list_name[time]:
            c+=1
        bodies[c].impact(data[time,0],data[time,1],data[time,2],data[time,3])

    for p in bodies:
        p.plot_thermal_history()

    plt.figure()
    plt.show()
