import numpy as np
import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
#import scienceplots
from filterpy.kalman import UnscentedKalmanFilter as UKF
from filterpy.kalman import MerweScaledSigmaPoints



def hovorka(t,y,par):
    G1,G2,Q1,Q2,S1,S2,I,x1,x2,x3,C=y

    (Vg, R_thr, R_cl, V1, ka, ke, Bio, tmax, kint,
    Ug_ceil, F01c, Fr, EGP, Ug, u, D, G,k12,EGP0) = par


# Interstitial glucose kinetics


    dCdt=kint*(G-C)
   

# glucose absorption

    dG1dt=(-G1/tmax)+(Bio*D)
    dG2dt=(G1/tmax)-(G2/tmax)


# glucose kinetics subsystem

    dQ1dt=-F01c-(x1*Q1)+(k12*Q2)-Fr+EGP0*(1-x3)+Ug
    dQ2dt=(x1*Q1)-(k12+x2)*Q2
    
#insulin absorption and kinetics subsystem

    dS1dt=u-(ka*S1)
    dS2dt=(ka*S1)-(ka*S2)
    dIdt=((ka*S2)/V1)-(ke*I)

# insulin action subsystem

    dx1dt=-(ka1*x1)+(S1t*kb1*I)
    dx2dt=-(ka2*x2)+(S1d*kb2*I)
    dx3dt=-(ka3*x3)+(S1e*kb3*I)


    

    return [dG1dt,dG2dt,dQ1dt,dQ2dt,dS1dt,dS2dt,dIdt,dx1dt,dx2dt,dx3dt,dCdt]


def next_state(state,par):
    t_span = (0,1) 
    #noise1=np.random.multivariate_normal(np.array([0,0,0]),Q)
    solution = solve_ivp(hovorka, t_span, state,args=(par,),method='RK45', max_step=0.1)

    y_value = solution.y[:,-1]

    #print(y_value)

    return y_value


def meal(n):
    
    if 60<=n<75:
        
        return 20
    elif 350<=n<375 :

        return 30
    #elif n==700 :
        #return 50

    #elif n==800:

        #return  50

    else:

        return 0
    


def hx(X):
    H=np.array([0,0,0,0,0,0,0,0,0,0,1])
    
    
    return np.array([np.dot(H.T,X)])




def fx(X,dt,par):
    t_span = (0,dt) 
    #noise1=np.random.multivariate_normal(np.array([0,0,0]),Q)
    solution = solve_ivp(hovorka, t_span,X,args=(par,),method='LSODA')

    y_value = solution.y[:,-1]

    #print(y_value)

    return y_value









#parameters

F01=.0097*70
F01s=F01/.85
EGP0=.016
k12=.066
S1d=.000505
S1e=.019
S1t=.001841
kb1=.005
kb2=.0005
kb3=.0023

ka1=.006
ka2=.06
ka3=.03


Vg=.16*70
R_thr=9
R_cl=.01
V1=.12*70
ka=.018
ke=.14
Bio=.8
tmax=40
temp=tmax
kint=1
Ug_ceil=.02




#initial
G_basal = 5
Q1 = G_basal * Vg
state=np.array([0,0,Q1,.5*Q1,0,0,0,.001,.001,.001,Q1*.1])



X0=state


dim_x=11
dim_z=1
dt=1

points = MerweScaledSigmaPoints(n=dim_x, alpha=0.1, beta=2, kappa=0)


ukf = UKF(dim_x=dim_x, dim_z=dim_z, dt=1,
          fx=lambda x, dt: fx(x, dt, par),
          hx=hx, points=points)

ukf.P = np.diag([1, 1, .5, .5, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 1])



ukf.Q = np.eye(dim_z)*.000001
ukf.R = np.array([[1]])






est1=[]
est2=[]

reading=[]





#u=0
#D=0
#EGP=16.9

ns=[]
ns1=[]

ns2=[]
ns3=[]
carb=[]
ins=[]


time=list(range(1000))



test=0
for i in range(len(time)):


    measurement=state[10]+np.random.normal(0,1)
    #measurement=(state[10])*(1-np.exp(-kint*i))+np.random.normal(0,.1)

    
   

    u=test
    
    #u=25

    #if 5<=i<=114:
        #u=200
         #u=2
        

    #elif 400<=i<=513:
        #u=200
         #u=2

    #elif 585<=i<590:
        #u=2
    #else:
       # u=500
    
    D= meal(i)
    

    carb.append(D)
    ins.append(u)
    
    
    Ug=state[1]/tmax
    G= max(state[2]/Vg,(35.69/Vg))
    
    ig=max(state[10]/Vg,.2)

    #print(G,ig)

    EGP=max((1-state[9])*EGP0,0)

    F01c=(F01s*G)/(G+1)
  
    if G>=R_thr:

        Fr=R_cl*(G-R_thr)*Vg
    else:
        Fr=0

    
    #tmax_ceil=state[1]/Ug_ceil


    #if Ug>Ug_ceil:

     #   tmax=tmax_ceil
    #else:
     #   tmax=temp

   

    

    par=(Vg, R_thr, R_cl, V1, ka, ke, Bio, tmax, kint,
    Ug_ceil, F01c, Fr, EGP, Ug, u, D, G,k12,EGP0)


    ukf.x[2] = max(ukf.x[2],35.69)
    #ukf.x[10]= max(ukf.x[10],3.69)
   
    

    est_g=ukf.x[2]/Vg
    est_ig=max(ukf.x[10]/Vg,.2)

    print(est_g)


    ukf.predict()

    ukf.update(measurement)
    
    #print(state[2])
    sol=next_state(state,par)

    #test=max(400-(400*np.exp(state[2]-sol[2])),0)
    test=max(sol[2]*np.exp(state[2]-sol[2]),0)+max(40*np.exp(-(state[9]-sol[9])),0)

    
    #print(ukf.x[2])

    est1.append(est_g)
    est2.append(est_ig)
   

    
    ns.append(sol[0])
    ns1.append(sol[1])
    ns2.append(G)
    ns3.append(ig)
    state=sol

    

    


    #reading.append(measurement/Vg)
    










#plt.style.use(['science', 'ieee'])

#plt.plot(time,ns,ns3)
#plt.show()


#plt.plot(time,ns1,ns2)
#plt.show()

#plt.plot(time,carb)
#plt.show()


fig, axs = plt.subplots(2,2, figsize=(6, 4))
axs[0,0].plot(time,ns2,label='truth',color='black', linestyle='solid')
axs[0,0].plot(time,est1,label='ukf',color='black', linestyle='dotted')
axs[0,0].legend()
axs[0,0].set_title('Blood glucose concentration ')
axs[0,0].set_ylabel('Glucose concentration (mmol/L)')
axs[0,0].set_xlabel('time (minutes)')



axs[0,1].plot(time,ns3,label='truth',color='black', linestyle='solid')
axs[0,1].plot(time,est2,label='ukf',color='black', linestyle='dotted')
axs[0,1].legend()
axs[0,1].set_title('Interstitial glucose concentration')
axs[0,1].set_ylabel('Glucose concentration (mmol/L)')
axs[0,1].set_xlabel('time (minutes)')
#axs[1].legend()

axs[1,0].plot(time,ins,color='black', linestyle='solid')
axs[1,0].set_title('insulin infusion rate')
axs[1,0].set_ylabel('insulin infusio rate (mU/minute)')
axs[1,0].set_xlabel('time (minutes)')

axs[1,1].plot(time,ns,label='gut',color='black', linestyle='solid')
axs[1,1].plot(time,ns1,label='intestine',color='black', linestyle='dotted')
axs[1,1].legend()
axs[1,1].set_title('Glucose mass in gut and intestine')
axs[1,1].set_ylabel('glucose mass(mmol)')
axs[1,1].set_xlabel('time (minutes)')






plt.tight_layout()



plt.show()



# normal bg range 3.9-5.5
# after meal <7.8


#minimize cost
#s.t
   #constraint 1   f1(S1)
   #constarint 2   f2(G)

#target= sum over full horizon 1/no of horizons[(x2*q2)*u + (1-u)(x1*q1)]^2
   # searching u in [2,25]

          #3.9<q1/vg<7.8

    







