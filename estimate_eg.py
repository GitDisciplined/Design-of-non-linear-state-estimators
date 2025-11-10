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
    Ug_ceil, F01c, Fr, EGP, Ug, u, D, G) = par
    
    # glucose absorption

    dG1dt=(-G1/tmax)+(Bio*D)
    dG2dt=(G1/tmax)-(G2/tmax)


# glucose kinetics subsystem

    dQ1dt=-F01c-(x1*Q1)+(k12*Q2)-Fr+EGP+Ug
    dQ2dt=(x1*Q1)-(k12+x2)*Q2
    
#insulin absorption and kinetics subsystem

    dS1dt=u-(ka*S1)
    dS2dt=(ka*S1)-(ka*S2)
    dIdt=((ka*S2)/V1)-(ke*I)

# insulin action subsystem

    dx1dt=-(ka1*x1)+(S1t*kb1*I)
    dx2dt=-(ka2*x2)+(S1d*kb2*I)
    dx3dt=-(ka3*x3)+(S1e*kb3*I)

# Interstitial glucose kinetics


    dCdt=kint*(G-C)
   
    

    return [dG1dt,dG2dt,dQ1dt,dQ2dt,dS1dt,dS2dt,dIdt,dx1dt,dx2dt,dx3dt,dCdt]


def next_state(state,par):
    t_span = (0,1) 
    #noise1=np.random.multivariate_normal(np.array([0,0,0]),Q)
    solution = solve_ivp(hovorka, t_span, state,args=(par,))

    y_value = solution.y[:,-1]

    #print(y_value)

    return y_value


def meal(n):
    
    if 5<=n<25:
        
        return 20
    #elif 200<=n<220  :

        #return 22
    #elif n==700 :
        #return 50

    #elif n==800:

        #return  50

    else:

        return 0
    


def hx(X):
    H=np.array([0,0,1,0,0,0,0,0,0,0,0])
    
    
    return np.array([np.dot(H.T,X)])




def fx(X,dt,par):
    t_span = (0,dt) 
    #noise1=np.random.multivariate_normal(np.array([0,0,0]),Q)
    solution = solve_ivp(hovorka, t_span,X,args=(par,))

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
kb1=.0034
kb2=.056
kb3=.024

ka1=S1t*kb1
ka2=S1d*kb2
ka3=S1e*kb3


Vg=.16*70
R_thr=9
R_cl=.01
V1=.12*70
ka=.018
ke=.14
Bio=.8
tmax=40
temp=tmax
kint=.07
Ug_ceil=.02




#initial
G_basal = 5
Q1 = G_basal * Vg
state=np.array([0,0,Q1,Q1*.1,0,0,0,0,0,0,5])
X0=state


dim_x=11
dim_z=1
dt=1

points = MerweScaledSigmaPoints(n=dim_x, alpha=0.001, beta=2, kappa=0)


ukf = UKF(dim_x=dim_x, dim_z=dim_z, dt=1.0,
          fx=lambda x, dt: fx(x, dt, par),
          hx=hx, points=points)

ukf.P = np.eye(dim_x)*.05
ukf.Q = np.eye(dim_x)*.0001
ukf.R = np.eye(dim_z)*10



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


time=list(range(1000))




for i in range(len(time)):


    #measure=state[10]+np.random.normal(5,1)
    measurement=(state[10])*(1-np.exp(-kint*i))+np.random.normal(0,.1)

    
   
    
    


    if 400<=i<=450:
        u=0
    elif 10<=i<=50:
        u=0
    else:
        u=0
    
    D= meal(i)
    

    carb.append(D)
    
    
    
    Ug=state[1]/tmax
    G= state[2]/Vg
    

    EGP=max((1-state[9])*EGP0,0)

    F01c=(F01s*G)/(G+1)
  
    if G>=R_thr:

        Fr=R_cl*(G-R_thr)*Vg
    else:
        Fr=0

    
    tmax_ceil=state[1]/Ug_ceil


    if Ug>Ug_ceil:

        tmax=tmax_ceil
    else:
        tmax=temp

   

    

    par=(Vg, R_thr, R_cl, V1, ka, ke, Bio, tmax, kint,
    Ug_ceil, F01c, Fr, EGP, Ug, u, D, G)


    ukf.predict()

    ukf.update(measurement)
    
    #print(state[2])
    sol=next_state(state,par)

    
    #print(ukf.x[2])

    est1.append(ukf.x[2]/Vg)
    est2.append(ukf.x[10])
   

    
    ns.append(sol[0])
    ns1.append(sol[1])
    ns2.append((sol[2]/Vg))
    ns3.append(sol[10])
    state=sol

    

    


    #reading.append(measurement/Vg)










#plt.style.use(['science', 'ieee'])

#plt.plot(time,ns,ns3)
#plt.show()


#plt.plot(time,ns1,ns2)
#plt.show()

#plt.plot(time,carb)
#plt.show()


fig, axs = plt.subplots(2, 2, figsize=(6, 4))
axs[0,0].plot(time,ns2)
axs[0,0].set_title('Blood glucose concentration ')
axs[0,0].set_ylabel('Glucose concentration (mmol/L)')
axs[0,0].set_xlabel('time (minutes)')



axs[0,1].plot(time,ns3,est2)
axs[0,1].set_title('Interstitial glucose concentration')
axs[0,1].set_ylabel('Glucose concentration (mmol/L)')
axs[0,1].set_xlabel('time (minutes)')
#axs[1].legend()

axs[1,0].plot(time,carb)
axs[1,0].set_title('meal intake rate')
axs[1,0].set_ylabel('Carb intake rate (mmol/minute)')
axs[1,0].set_xlabel('time (minutes)')

axs[1,1].plot(time,ns,ns1)
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

    







