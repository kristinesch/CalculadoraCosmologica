import numpy as np
import scipy as sc
import scipy.integrate as integrate
import scipy.constants as const
import math
from matplotlib import pyplot as plt

#TODO fiks dL for neg rot
#Funker jo ikke for omega_K lik 0 ??

def E(z,Om_R, Om_M, Om_Lambda, Om_K): #z is a np array
    pM=np.power((1+z),3)
    pR=np.power((1+z),4)
    pK=np.power((1+z),2)
    return 1/(np.sqrt(Om_M*pM+Om_R*pR+Om_K*pK+Om_Lambda))

def integrateE(z,Om_R, Om_M, Om_Lambda, Om_K):
    Eintegrated=[]
    for zi in z:
        x=integrate.quad(E,z[0],zi,(Om_R, Om_M, Om_Lambda, Om_K))
        Eintegrated.append(x[0])
    return np.array(Eintegrated)


def calculate_dL(z,Om_R, Om_M, Om_Lambda, Om_K,H0):
    #Eintegrated=integrate.quad(E,z[0],z[-1],(Om_R, Om_M, Om_Lambda, Om_K))
    Eintegrated=integrateE(z,Om_R, Om_M, Om_Lambda, Om_K)
    print("E:",Eintegrated)
    #print(const.c*(1+z)/(H0*math.sqrt(Om_K)))
    print(np.sinh(Eintegrated))
    if Om_K>=0:
        return const.c*(1+z)/(H0*math.sqrt(Om_K))*np.sinh(math.sqrt(Om_K)*Eintegrated)
    else: #squareroot gives imaginary number, with i counteracted by the fact that sinh(ix)=isin(x)
        #return math.sqrt(abs(Om_K))*Eintegrated
        return const.c*(1+z)/(H0*math.sqrt(abs(Om_K)))*np.sin(math.sqrt(abs(Om_K))*Eintegrated) #denne er det noe rart med... skal vel egt være minustegn?


def calculate_dA(z,dL):
    return (1+np.power(z,-2))*dL

def plot_distance(z,ds,filename,dstring,Om_M,Om_R,Om_Lambda,H0): #ds a list
    #colors=[""]
    fig,ax=plt.subplots(1,1)
    for i, di in enumerate(ds):
        ax.plot(z,di,label="$\Omega_M=$,"+str(Om_M[i])+" $\Omega_R=$,"+str(Om_R[i])+"$ \Omega_\Lambda=$,"+str(Om_Lambda[i])+" $H_0$="+str(H0[i]))
    ax.set_xlabel("z")
    ax.set_ylabel(dstring)
    #plt.show()
    fig.legend()
    fig.savefig(filename)
    plt.show()

def distances(z,Om_R, Om_M, Om_Lambda, H0,k):
    #Om_K=(-const.c*const.c*k)/(H0*H0)
    Om_K=1-Om_Lambda-Om_M-Om_R
    print("K:",Om_K)

    dL=calculate_dL(z,Om_R, Om_M, Om_Lambda, Om_K,H0)
    dA=calculate_dA(z,dL)

    #plot_distance(z,dL,"dLfig","dL")
    #plot_distance(z,dA,"dAfig","dA")
    return dL,dA

def main():

    #Take as input?
    H0=[67,67,67,67,67]
    Om_M=[0,0.3,0.5,0,0]
    Om_R=[0,0,0,0.5,0]
    Om_Lambda=[0,0,0,0,0.5]

    k=1

    N=10
    z_min=0
    z_max=15
    z=np.linspace(z_min,z_max,N)

    #running distances etc
    dLs=[]
    dAs=[]   
    for i in range(len(H0)):
        print(i)
        dLi,dAi= distances(z,Om_R[i], Om_M[i], Om_Lambda[i], H0[i],k)
        dLs.append(dLi)
        dAs.append(dAi)

    plot_distance(z,dAs,"dAfig","dA",Om_M,Om_R,Om_Lambda,H0)
    plot_distance(z,dLs,"dLfig","dL",Om_M,Om_R,Om_Lambda,H0)



main()

