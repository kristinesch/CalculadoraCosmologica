import numpy as np
import scipy as sc
import scipy.integrate as integrate
import scipy.constants as const
import math
from matplotlib import pyplot as plt

#TODO:
#remove k
#fix bug for odeint
#fix H in legend

""" a(t) and H(t)-----------------------------------------------------------------"""
def g(y, t, H0, Om_M, Om_R, Om_Lambda):
    a, aprime = y
    dydt = [aprime, -(a/2)*H0**2*(Om_M/(a**3)+2*Om_R/(a**4)-2*Om_Lambda)]
    return dydt


def calculate_a(y0,t1,t2,H0,Om_M,Om_R,Om_Lambda):

    sol1 = integrate.odeint(g, y0, t1, args=(H0, Om_M, Om_R, Om_Lambda))
    sol2 = integrate.odeint(g, y0, t2, args=(H0, Om_M, Om_R, Om_Lambda))
    return sol1,sol2

def g2(t,a,c):
    return c/a



def calculate_dP(t,a,c):
    aint=[]
    for ti in t:
        x=integrate.quad(g2,0,ti,args=(a,c))
        aint.append(x[0])
    return a*np.array(aint)


def time_params(y0,t1,t2,H0List,Om_MList,Om_RList,Om_LambdaList,Hlist):
    c=const.c
    sol1List=[]
    sol2List=[]
    colors=["red","blue","green","orange","purple","black"]
    for i in range(len(H0List)):
        sol1,sol2=calculate_a(y0,t1,t2,H0List[i],Om_MList[i],Om_RList[i],Om_LambdaList[i])
        sol1List.append(sol1)
        sol2List.append(sol2)

    figa, axa=plt.subplots(1,1)
    figH, axH=plt.subplots(1,1)
    figdH, axdH=plt.subplots(1,1)
    figdP, axdP=plt.subplots(1,1)

    for i in range(len(sol1List)):
        a1=sol1List[i][:,0]
        a2=sol2List[i][:,0]
        adot1=sol1List[i][:,1]
        adot2=sol2List[i][:,1]

        #x = x[~numpy.isnan(x)] remove nans?

        H1=(adot1/a1)*(3.086e22/(1e3*3600*24*365*10**9))
        H2=(adot2/a2)*(3.086e22/(1e3*3600*24*365*10**9))
        Mpc=3.0857e22
        dH1=c/(H1)
        dH2=c/(H2)

        #dH1 = dH1[~np.isnan(dH1)]
        #dH2 = dH2[~np.isnan(dH2)]

        dP1=calculate_dP(t1,a1,c)
        #dP2=calculate_dP(t2,a2,c)

        axa.plot(t1,a1, color=colors[i],label="$\Omega_M=$"+str(Om_MList[i])+" $, \Omega_R=$"+str(Om_RList[i])+
        "$, \Omega_\Lambda=$"+str(Om_LambdaList[i])+" $, H_0$="+str(Hlist[i]))
        axa.plot(t2,a2,color=colors[i])

        axH.plot(t1,H1, color=colors[i],label="$\Omega_M=$"+str(Om_MList[i])+" $, \Omega_R=$"+str(Om_RList[i])+
        "$, \Omega_\Lambda=$"+str(Om_LambdaList[i])+" $, H_0$="+str(Hlist[i]))
        axH.plot(t2,H2,color=colors[i])

        axdH.plot(t1,dH1, color=colors[i],label="$\Omega_M=$"+str(Om_MList[i])+" $, \Omega_R=$"+str(Om_RList[i])+
        "$, \Omega_\Lambda=$"+str(Om_LambdaList[i])+" $, H_0$="+str(Hlist[i]))
        axdH.plot(t2,dH2,color=colors[i])

        axdP.plot(t1,dP1, color=colors[i],label="$\Omega_M=$"+str(Om_MList[i])+" $, \Omega_R=$"+str(Om_RList[i])+
        "$, \Omega_\Lambda=$"+str(Om_LambdaList[i])+" $, H_0$="+str(Hlist[i]))
        axdH.plot(t2,dH2,color=colors[i])

    figa.suptitle("a(t)")
    axa.set_xlabel("t [Gyr]")
    axa.set_ylabel("a(t)")
    axa.set_ylim(0.00001,6)
    figa.legend(loc="center")
    figa.savefig("a_fig")
    

    figH.suptitle("H(t)")
    axH.set_xlabel("t [Gyr]")
    axH.set_ylabel("H(t) $[kms^{-1}Mpc^{-1}]$")
    axH.set_ylim(0,300)
    figH.legend(loc="center")
    figH.savefig("H_fig")

    figdH.suptitle("$d_H(t)$")
    axdH.set_xlabel("t [Gyr]")
    axdH.set_ylabel("$d_H(t)$")
    axdH.set_yscale("log")
    axdH.set_ylim(1e4,1e09)
    figdH.legend(loc="lower right")
    figdH.savefig("dH_fig")

    figdP.suptitle("$d_P(t)$")
    axdP.set_xlabel("t [Gyr]")
    axdP.set_ylabel("$d_P(t)$")
    axdP.set_yscale("log")
    figdP.legend(loc="lower right")
    figdP.savefig("dP_fig")

    plt.show()




"""dL and dA------------------------------------------------------------------"""

def E(z,Om_R, Om_M, Om_Lambda, Om_K): 
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
    Eintegrated=integrateE(z,Om_R, Om_M, Om_Lambda, Om_K)
    if Om_K>0:
        return const.c*(1+z)/(H0*math.sqrt(Om_K))*np.sinh(math.sqrt(Om_K)*Eintegrated)
    elif Om_K==0:
        return const.c*(1+z)/H0*Eintegrated
    else: #squareroot gives imaginary number, with i counteracted by the fact that sinh(ix)=isin(x)
        print("neg")
        return const.c*(1+z)/(H0*math.sqrt(abs(Om_K)))*np.sin(math.sqrt(abs(Om_K))*Eintegrated) 


def calculate_dA(z,dL):
    return (np.power(1+z,-2))*dL

def plot_distance(z,ds,filename,dstring,Om_M,Om_R,Om_Lambda,Hlist): 
    colors=["red","blue","green","orange","purple","black"]
    fig,ax=plt.subplots(1,1)
    for i, di in enumerate(ds):
        ax.plot(z,di,color=colors[i],label="$\Omega_M=$"+str(Om_M[i])+" $,\Omega_R=$"+str(Om_R[i])+"$, \Omega_\Lambda=$"+str(Om_Lambda[i])+" $, H_0$="+str(Hlist[i]))
    ax.set_xlabel("z")
    ax.set_ylabel(dstring)
    fig.legend()
    fig.savefig(filename)
    plt.show()

def calculate_dA_dL(z,Om_R, Om_M, Om_Lambda, H0,k):
    Om_K=1-Om_Lambda-Om_M-Om_R
    dL=calculate_dL(z,Om_R, Om_M, Om_Lambda, Om_K,H0)
    dA=calculate_dA(z,dL)
    return dL,dA


def redshift_params(z,Om_R, Om_M, Om_Lambda, H0,k,Hlist):
    dLs=[]
    dAs=[]   
    for i in range(len(H0)):
        dLi,dAi= calculate_dA_dL(z,Om_R[i], Om_M[i], Om_Lambda[i], H0[i],k)
        dLs.append(dLi)
        dAs.append(dAi)

    plot_distance(z,dAs,"dAfig","dA",Om_M,Om_R,Om_Lambda,Hlist)
    plot_distance(z,dLs,"dLfig","dL",Om_M,Om_R,Om_Lambda,Hlist)


"""--------------------------------------------------------------------"""



def main():
    year=365*24*60*60
    h=0.7
    H=h*100
    H0=H/(3.0857*1e19)*1e9*year #unit: per Gy
    H0s=[H0,H0,H0,H0,H0,0.8*H0]
    Om_Ms=[0.3,0.3,5,1,0,0.3]
    Om_Rs=[0,0,0,0,1,0]
    Om_Lambdas=[0.7,0,0,0,0,0.7]
    y0=np.array([1,H0])
    
    k=0
    Hlist=[H,H,H,H,H,0.8*H]

    Nz=1000
    z_min=0
    z_max=15
    z=np.linspace(z_min,z_max,Nz)

    Nt=1000
    t_min=-14.15
    t_max=30
    t1 = np.linspace(0, t_max, Nt)
    t2 = np.linspace(0, t_min, Nt)


    redshift_params(z,Om_Rs, Om_Ms, Om_Lambdas, H0s,k,Hlist)


    time_params(y0,t1, t2, H0s, Om_Ms,Om_Rs,Om_Lambdas,Hlist)





main()


