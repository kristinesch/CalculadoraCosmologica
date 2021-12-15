import numpy as np
import scipy as sc
import scipy.integrate as integrate
import scipy.constants as const
import math
from matplotlib import pyplot as plt

#TODO fiks dL for neg rot
#Funker jo ikke for omega_K lik 0 ??

""" a(t) and H(t)-----------------------------------------------------------------"""

def f(x,t,H0,Om_M,Om_R,Om_Lambda):
    a=x[0]
    #return np.array([x[1],-(a/2)*H0*H0*(Om_M/(np.power(a,3))+2*Om_R/(np.power(a,4))-2*Om_Lambda)])
    return np.array([x[1],-(1/2)*H0*H0*(Om_M/(np.power(a,1))+2*Om_R/(np.power(a,2))-2*Om_Lambda)])
#testing....
def g(x,t,b,c):
    a=x[0]
    return np.array([x[1],-a*c*(b/np.power(a,3)+c)])
    
def calculate_a_H(t,H0,Om_M,Om_R,Om_Lambda):
    x0=np.array([H0,1])
    x=integrate.odeint(f,x0,t, args=(H0,Om_M,Om_R,Om_Lambda))
    #b=2
    #c=2
    #x=integrate.odeint(g,x0,t,args=(b,c))
    #print(x)
    a=x[:,0]
    adot=x[:,1]
    H=a/adot
    return a, H

def plot_a(t,aList,H0s,Om_Ms,Om_Rs,Om_Lambdas,filename):
    fig, ax=plt.subplots(1,1)
    for i, ai in enumerate(aList):
        ax.plot(t,ai, label="$\Omega_M=$,"+str(Om_Ms[i])+" $\Omega_R=$,"+str(Om_Rs[i])+"$ \Omega_\Lambda=$,"+str(Om_Lambdas[i])+" $H_0$="+str(H0s[i]))
    ax.set_xlabel("t")
    ax.set_ylabel("a")
    #plt.show()
    fig.legend()
    fig.savefig(filename)
    plt.show() 

def plot_H(t,HList,H0s,Om_Ms,Om_Rs,Om_Lambdas,filename):
    fig, ax=plt.subplots(1,1)
    for i, Hi in enumerate(HList):
        ax.plot(t,Hi, label="$\Omega_M=$,"+str(Om_Ms[i])+" $\Omega_R=$,"+str(Om_Rs[i])+"$ \Omega_\Lambda=$,"+str(Om_Lambdas[i])+" $H_0$="+str(H0s[i]))
    ax.set_xlabel("t")
    ax.set_ylabel("H")
    #plt.show()
    fig.legend()
    fig.savefig(filename)
    plt.show() 

def a_and_H(t,H0s,Om_Ms,Om_Rs,Om_Lambdas):
    HList=[]
    aList=[]
    for i in range(len(H0s)):
        a,H=calculate_a_H(t,H0s[i],Om_Ms[i],Om_Rs[i],Om_Lambdas[i])
        aList.append(a)
        HList.append(H)
    plot_a(t,aList,H0s,Om_Ms,Om_Rs,Om_Lambdas,"a_fig")
    plot_H(t,HList,H0s,Om_Ms,Om_Rs,Om_Lambdas,"H_fig")


def test():
    t=np.linspace(0,10,1000)
    H0=67
    Om_M=1
    Om_R=0
    Om_Lambda=0
    a,H=calculate_a_H(t,H0,Om_M,Om_R,Om_Lambda)
    print(a)
    plt.plot(t,H)
    plt.show()
    plt.plot(t,a)
    plt.show()

#test()



"""dL and dA------------------------------------------------------------------"""

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

def calculate_dA_dL(z,Om_R, Om_M, Om_Lambda, H0,k):
    #Om_K=(-const.c*const.c*k)/(H0*H0)
    Om_K=1-Om_Lambda-Om_M-Om_R
    print("K:",Om_K)

    dL=calculate_dL(z,Om_R, Om_M, Om_Lambda, Om_K,H0)
    dA=calculate_dA(z,dL)

    #plot_distance(z,dL,"dLfig","dL")
    #plot_distance(z,dA,"dAfig","dA")
    return dL,dA


def distances(z,Om_R, Om_M, Om_Lambda, H0,k):
    #running distances etc
    dLs=[]
    dAs=[]   
    for i in range(len(H0)):
        dLi,dAi= calculate_dA_dL(z,Om_R[i], Om_M[i], Om_Lambda[i], H0[i],k)
        dLs.append(dLi)
        dAs.append(dAi)

    plot_distance(z,dAs,"dAfig","dA",Om_M,Om_R,Om_Lambda,H0)
    plot_distance(z,dLs,"dLfig","dL",Om_M,Om_R,Om_Lambda,H0)


"""--------------------------------------------------------------------"""



def main():

    #Take as input?
    H0s=[67,67,67,67,67]
    Om_Ms=[0,0.3,0.5,0,0.4]
    Om_Rs=[0,0,0,0.5,0]
    Om_Lambdas=[0,0,0,0,0.5]

    k=1

    Nz=10
    z_min=0
    z_max=15
    z=np.linspace(z_min,z_max,Nz)

    Nt=10000
    year=365*24*60*60
    t_min=0
    t_max=1000
    t=np.linspace(t_min,t_max,Nt)

    #distances(z,Om_Rs, Om_Ms, Om_Lambdas, H0s,k)
    a_and_H(t, H0s, Om_Ms,Om_Rs,Om_Lambdas)





main()



# def k(x,t,a,b):
#     return np.array([x[1],a+b])

# t=np.linspace(0,10,100)
# a=2
# b=0
# sol=integrate.odeint(k,[0,0],t,args=(a,b))
# plt.plot(t,sol[:,0])
# plt.show()