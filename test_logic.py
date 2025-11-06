import numpy as np
import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd
from filterpy.kalman import UnscentedKalmanFilter as UKF
from filterpy.kalman import MerweScaledSigmaPoints
#import scienceplots



#parameters

F01=11.1
F01s=F01/.85 
EGP0=.016   
k12=.06
S1d=.000505
S1e=.019
S1t=.001841
kb1=.0034
kb2=.056
kb3=.024

ka1=S1t*kb1
ka2=S1d*kb2
ka3=S1e*kb3


Vg=10
R_thr=9
R_cl=.01
V1=8
ka=.018
ke=.14
Bio=.8
tmax=40
temp=tmax
kint=.05
Ug_ceil=.02



def hovorka(t,y,par):

    G1,G2,Q1,Q2,S1,S2,I,x1,x2,x3=y

    (Vg, R_thr, R_cl, V1, ka, ke, Bio, tmax,
    Ug_ceil, F01c, Fr, EGP, Ug, u, D, G,kint,
     F01,F01s,S1d,S1e,S1t,kb1,kb2,kb3,ka1,ka2,ka3) = par

    
    
    # glucose absorption

    dG1dt=(-G1/tmax)+((Bio*D))
    dG2dt=(G1/tmax)-Ug


# glucose kinetics subsystem

    dQ1dt=-F01c-(x1*Q1)+(k12*Q2)-Fr+EGP+Ug
    dQ2dt=(x1*Q1)-(k12+x2)*Q2
    
#insulin absorption and kinetics subsystem

    dS1dt=u-(ka*S1)
    dS2dt=(ka*S1)-(ka*S2)
    dIdt=((ka*S2)/V1)-(ke*I)

# insulin action subsystem

    dx1dt=-(ka1*x1)+(kb1*I)
    dx2dt=-(ka2*x2)+(kb2*I)
    dx3dt=-(ka3*x3)+(kb3*I)

# Interstitial glucose kinetics


    #dCdt=kint*(G-C)
   
    

    return [dG1dt,dG2dt,dQ1dt,dQ2dt,dS1dt,dS2dt,dIdt,dx1dt,dx2dt,dx3dt]


def fx(X,dt,par):
    t_span = (0,dt) 
    #noise1=np.random.multivariate_normal(np.array([0,0,0]),Q)
    solution = solve_ivp(hovorka, t_span,X,args=(par,))

    y_value = solution.y[:,-1]

    #print(y_value)

    return y_value



def meal(n):
    
    if 185<=n<200 :

        return 25
        

    else:

        return 0
    
def hx(X):
    H=np.array([0,0,1,0,0,0,0,0,0,0])
    
    
    return np.array([np.dot(H.T,X)])



dim_x=10
dim_z=1
dt=10


#df = pd.read_csv("results.csv")
#ig = df['CGM'].to_numpy()


X0=np.array([0,0,50,0,300,0,0,0,0,0])

points = MerweScaledSigmaPoints(n=dim_x, alpha=0.001, beta=2, kappa=7)

ukf = UKF(dim_x=dim_x, dim_z=dim_z, dt=1.0,
          fx=lambda x, dt: fx(x, dt, par),
          hx=hx, points=points)

ukf.P = np.eye(dim_x)*1
ukf.Q = np.eye(dim_x)*1
ukf.R = np.eye(dim_z) *.01
#steps=len(ig)

steps=1000


estimate=[]
estimate2=[]
reading=[]

time=np.linspace(0,steps-1,steps)


ukf.x=X0

for  itere in range(steps):
    D=meal(itere)
    u=0
    measurement=ukf.x[2]

    Ug=ukf.x[1]/tmax
    G= ukf.x[2]/Vg
    

    EGP=max((1+ukf.x[9])*EGP0,0)

    F01c=(F01s*G)/(G+1)
  
    if G>=R_thr:

        Fr=R_cl*(G-R_thr)*Vg
    else:
        Fr=0

    
    tmax_ceil=ukf.x[1]/Ug_ceil


    if Ug>Ug_ceil:

        tmax=tmax_ceil
    else:
        tmax=temp





    #if itere==90:
       # u=.000006
    #else:
        #u=0

    #if G>=5:
        
        #u=25
    #else:
        #u=0
  
        

    par=(Vg, R_thr, R_cl, V1, ka, ke, Bio, tmax,
    Ug_ceil, F01c, Fr,EGP, Ug, u, D, G,kint,
     F01,F01s,S1d,S1e,S1t,kb1,kb2,kb3,ka1,ka2,ka3)

    #ukf.x = X0

    ukf.predict()
    #print(ukf.x)

    #measurement=ig[itere]+np.random.normal(0,4)
   
    ukf.update(measurement)

    #X0=ukf.x

    #print(ukf.x)

    estimate.append((ukf.x[2]/Vg))
    #estimate2.append(ukf.x[0])
    reading.append((measurement/Vg))


plt.plot(time,estimate,reading)
plt.show()








    



   
 





