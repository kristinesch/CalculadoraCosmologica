import numpy as np
import scipy as sc
import scipy.integrate as integrate
import scipy.constants as const
import math
from matplotlib import pyplot as plt

#TODO fiks dL for neg rot
#Funker jo ikke for omega_K lik 0 ??
#Fjerne k

""" a(t) and H(t)-----------------------------------------------------------------"""
def g(y, t, H0, Om_M, Om_R, Om_Lambda):
    a, aprime = y
    dydt = [aprime, -(a/2)*H0**2*(Om_M/(a**3)+2*Om_R/(a**4)-2*Om_Lambda)]
    return dydt

# def f(x,t,H0,Om_M,Om_R,Om_Lambda):
#     return np.array([x[1],-(x[0]/2)*H0*H0*(Om_M/(np.power(x[0],3))+2*Om_R/(np.power(x[0],4))-2*Om_Lambda)])

def calculate_a(y0,t1,t2,H0,Om_M,Om_R,Om_Lambda):

    sol1 = integrate.odeint(g, y0, t1, args=(H0, Om_M, Om_R, Om_Lambda))
    sol2 = integrate.odeint(g, y0, t2, args=(H0, Om_M, Om_R, Om_Lambda))

    # a=np.append(sol2[:,0],sol1[:,0])
    # adot=np.append(sol2[:,1],sol1[:,1])
    # t=np.append(t2,t1)
    # print(t,a,adot)

    # plt.plot(t1, sol1[:, 0], 'b', t2, sol2[:,0], 'r', label='a(t)')
    # #plt.plot(t,a, 'r', label='a(t)')
    # # plt.plot(t, sol[:, 1], 'g', label='omega(t)')
    # plt.legend(loc='best')
    # plt.xlabel('t')
    # plt.grid()
    # plt.show()
    return sol1,sol2

def a_and_H(y0,t1,t2,H0List,Om_MList,Om_RList,Om_LambdaList):
    c=const.c
    sol1List=[]
    sol2List=[]
    colors=["red","blue","green","orange"]
    for i in range(len(H0List)):
        sol1,sol2=calculate_a(y0,t1,t2,H0List[i],Om_MList[i],Om_RList[i],Om_LambdaList[i])
        sol1List.append(sol1)
        sol2List.append(sol2)

    figa, axa=plt.subplots(1,1)
    figH, axH=plt.subplots(1,1)
    figdH, axdH=plt.subplots(1,1)

    for i in range(len(sol1List)):
        a1=sol1List[i][:,0]
        a2=sol2List[i][:,0]
        adot1=sol1List[i][:,1]
        adot2=sol2List[i][:,1]

        #x = x[~numpy.isnan(x)] remove nans?

        H1=(adot1/a1)*(3.086e22/(1e3*3600*24*365*10**9))
        H2=(adot2/a2)*(3.086e22/(1e3*3600*24*365*10**9))
        dH1=c/H1
        dH2=c/H2

        #dH1 = dH1[~np.isnan(dH1)]
        #dH2 = dH2[~np.isnan(dH2)]

        axa.plot(t1,a1, color=colors[i],label="$\Omega_M=$"+str(Om_MList[i])+" $, \Omega_R=$"+str(Om_RList[i])+
        "$, \Omega_\Lambda=$"+str(Om_LambdaList[i])+" $, H_0$="+str(round(H0List[i],2)))
        axa.plot(t2,a2,color=colors[i])

        axH.plot(t1,H1, color=colors[i],label="$\Omega_M=$"+str(Om_MList[i])+" $, \Omega_R=$"+str(Om_RList[i])+
        "$, \Omega_\Lambda=$"+str(Om_LambdaList[i])+" $, H_0$="+str(round(H0List[i],2)))
        axH.plot(t2,H2,color=colors[i])

        axdH.plot(t1,dH1, color=colors[i],label="$\Omega_M=$"+str(Om_MList[i])+" $, \Omega_R=$"+str(Om_RList[i])+
        "$, \Omega_\Lambda=$"+str(Om_LambdaList[i])+" $, H_0$="+str(round(H0List[i],2)))
        axdH.plot(t2,dH2,color=colors[i])

    figa.suptitle("a(t)")
    axa.set_xlabel("t [Gyr]")
    axa.set_ylabel("a(t)")
    axa.set_ylim(0.00001,6)
    figa.legend(loc="center")
    figa.savefig("a_fig")
    

    figH.suptitle("H(t)")
    axH.set_xlabel("t [Gyr]")
    axH.set_ylabel("H(t)")
    axH.set_ylim(0,300)
    figH.legend(loc="center")
    figH.savefig("H_fig")

    figdH.suptitle("$d_H(t)$")
    axdH.set_xlabel("t [Gyr]")
    axdH.set_ylabel("$d_H(t)$")
    axdH.set_yscale("log")
    axdH.set_ylim(1e04,1e09)
    figdH.legend(loc="lower right")
    figdH.savefig("dH_fig")

    plt.show()


"""a and H failed attempt"""

# def f(x,t,H0,Om_M,Om_R,Om_Lambda):
#     #a=x[0]
#     return np.array([x[1],-(x[0]/2)*H0*H0*(Om_M/(np.power(x[0],3))+2*Om_R/(np.power(x[0],4))-2*Om_Lambda)])

# def g(t,x,H0,Om_M,Om_R,Om_Lambda):
#     return np.array([x[1],-(x[0]/2)*H0*H0*(Om_M/(np.power(x[0],3)+10e-10)+2*Om_R/(np.power(x[0],4)+10e-10)-2*Om_Lambda)])
    
# def calculate_a_H(t,H0,Om_M,Om_R,Om_Lambda):
#     x0=np.array([1,H0])
#     #x=integrate.odeint(f,x0,t, args=(H0,Om_M,Om_R,Om_Lambda))
#     sol=integrate.solve_ivp(g,[t[0],t[-1]],x0, args=(H0,Om_M,Om_R,Om_Lambda),method="LSODA")
#     #print(sol)
#     t=sol.t
#     a=sol.y[0]
#     adot=sol.y[1]
#     print(t,a)
#     # a=x[:,0]
#     # adot=x[:,1]
#     H=a/adot
    
#     return t,a, H

# def plot_a(tList1,tList2,aList1,aList2,H0s,Om_Ms,Om_Rs,Om_Lambdas,filename):
#     fig, ax=plt.subplots(1,1)
#     for i, ai in enumerate(aList1):
#         ax.plot(tList1[i],ai, label="$\Omega_M=$"+str(Om_Ms[i])+" $, \Omega_R=$"+str(Om_Rs[i])+"$, \Omega_\Lambda=$"+str(Om_Lambdas[i])+" $, H_0$="+str(H0s[i]))
#         ax.plot(tList2[i],aList2[i])
#     ax.set_xlabel("t")
#     ax.set_ylabel("a")
#     fig.legend()
#     fig.savefig(filename)
#     plt.show() 

# def plot_H(tList1,tList2,HList1,HList2,H0s,Om_Ms,Om_Rs,Om_Lambdas,filename):
#     fig, ax=plt.subplots(1,1)
#     for i, Hi in enumerate(HList1):
#         ax.plot(tList1[i],Hi, label="$\Omega_M=$"+str(Om_Ms[i])+" $, \Omega_R=$"+str(Om_Rs[i])+"$, \Omega_\Lambda=$"+str(Om_Lambdas[i])+" $, H_0$="+str(H0s[i]))
#         ax.plot(tList2,HList2[i])
#     ax.set_xlabel("t")
#     ax.set_ylabel("H")
#     fig.legend()
#     fig.savefig(filename)
#     plt.show() 

# def a_and_H(t1,t2,H0s,Om_Ms,Om_Rs,Om_Lambdas):
#     HList1=[]
#     aList1=[]
#     tList1=[]

#     HList2=[]
#     aList2=[]
#     tList2=[]
#     for i in range(len(H0s)):
#         t,a,H=calculate_a_H(t1,H0s[i],Om_Ms[i],Om_Rs[i],Om_Lambdas[i])
#         aList1.append(a)
#         HList1.append(H)
#         tList1.append(t)

#         t,a,H=calculate_a_H(t2,H0s[i],Om_Ms[i],Om_Rs[i],Om_Lambdas[i])
#         aList2.append(a)
#         HList2.append(H)
#         tList2.append(t)


#     plot_a(tList1,tList2,aList1,aList2,H0s,Om_Ms,Om_Rs,Om_Lambdas,"a_fig")
#     plot_H(tList1,tList2,HList1,HList2,H0s,Om_Ms,Om_Rs,Om_Lambdas,"H_fig")




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
    #print("E:",Eintegrated)
    #print(const.c*(1+z)/(H0*math.sqrt(Om_K)))
    #print(np.sinh(Eintegrated))
    if Om_K>0:
        print("pos")
        return const.c*(1+z)/(H0*math.sqrt(Om_K))*np.sinh(math.sqrt(Om_K)*Eintegrated)
    elif Om_K==0:
        print("null")
        return const.c*(1+z)/H0*Eintegrated
    else: #squareroot gives imaginary number, with i counteracted by the fact that sinh(ix)=isin(x)
        #return math.sqrt(abs(Om_K))*Eintegrated
        print("neg")
        return const.c*(1+z)/(H0*math.sqrt(abs(Om_K)))*np.sin(math.sqrt(abs(Om_K))*Eintegrated) #denne er det noe rart med... skal vel egt være minustegn?


def calculate_dA(z,dL):
    return (np.power(1+z,-2))*dL

def plot_distance(z,ds,filename,dstring,Om_M,Om_R,Om_Lambda,H0): #ds a list
    #colors=[""]
    fig,ax=plt.subplots(1,1)
    for i, di in enumerate(ds):
        ax.plot(z,di,label="$\Omega_M=$"+str(Om_M[i])+" $,\Omega_R=$"+str(Om_R[i])+"$, \Omega_\Lambda=$"+str(Om_Lambda[i])+" $, H_0$="+str(H0[i]))
    ax.set_xlabel("z")
    ax.set_ylabel(dstring)
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
    year=365*24*60*60
    #Take as input?
    H0=1000/3.086e22*(3600*24*365*10**9)*67
    
    #H0=67
    H0s=[H0,H0,H0,H0]
    Om_Ms=[0.3,0.3,5,1]
    Om_Rs=[0,0,0,0]
    Om_Lambdas=[0,0.7,0,0]

    k=0

    Nz=1000
    z_min=0
    z_max=15
    z=np.linspace(z_min,z_max,Nz)

    Nt=10000
    
    t_min=-14.15
    t_max=30
    t1 = np.linspace(0, t_max, Nt)
    t2 = np.linspace(0, t_min, Nt)

    y0=np.array([1,H0])

    distances(z,Om_Rs, Om_Ms, Om_Lambdas, H0s,k)



    #a_and_H(y0,t1, t2, H0s, Om_Ms,Om_Rs,Om_Lambdas)





main()


#testing odeint...
# def k(x,t,a,b):
#     return np.array([x[1],a+b])

# t=np.linspace(0,10,100)
# a=2
# b=0
# sol=integrate.odeint(k,[0,0],t,args=(a,b))
# plt.plot(t,sol[:,0])
# plt.show()

a=np.array([1,2])
b=np.array([3,4])
print(np.append(a,b))

print(1000/3.086e22*(3600*24*365*10**(9))*67)