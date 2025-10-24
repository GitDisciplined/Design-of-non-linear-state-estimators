import numpy as np
import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def hovorka(t,y):
    G1,G2,Q1,Q2,S1,S2,I,x1,x2,x3,C=y

    
# glucose absorption

    dG1dt=(-G1/tmax)+(Bio*D)
    dG2dt=(G1/tmax)-(G2/tmax)


# glucose kinetics subsystem

    dQ1dt=-(F01c/(Vg*G))*Q1+(k12*Q2)-Fr+EGP+Ug
    dQ2dt=(x1*Q1)-(k12+x2)*Q2
    
#insulin absorption and kinetics subsystem

    dS1dt=u-(ka*S1)
    dS2dt=(ka*S1)-(ka*S2)
    dIdt=((ka*S2)/V1)-ke*I

# insulin action subsystem

    dx1dt=-ka1*x1+(S1t*kb1*I)
    dx2dt=-ka2*x2+(S1d*kb2*I)
    dx3dt=-ka3*x3+(S1e*kb3*I)

# Interstitial glucose kinetics


    dCdt=kint*(G-C)
   
    

    return [dG1dt,dG2dt,dQ1dt,dQ2dt,dS1dt,dS2dt,dIdt,dx1dt,dx2dt,dx3dt,dCdt]


def next_state(state):
    t_span = (0, 1) 
    #noise1=np.random.multivariate_normal(np.array([0,0,0]),Q)
    solution = solve_ivp(hovorka, t_span, state)

    y_value = solution.y[:,-1]

    #print(y_value)

    return y_value



#parameters

F01=11.1
F01s=F01/.85
EGP0=16.9
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


Vg=.16
R_thr=9
R_cl=.01
V1=.12
ka=.018
ke=.14
Bio=70
tmax=1/.3
kint=.066
Ug_ceil=.02



#initial

state=np.array([0,0,140,160,0,0,0,0,0,0,0])
u=0
D=5
EGP=0
ns=[]
time=list(range(200))

for i in range(len(time)):
    

    G= state[2]/Vg

    if EGP>=0:
        EGP=(1+state[9])*EGP0
    else:
        EGP=0
        

    F01c=(F01s*G)/(G+1)

    if G>=R_thr:

        Fr=R_cl*(G-R_thr)*Vg
    else:
        Fr=0

    Ug=state[1]/tmax

    tmax_ceil=state[1]/Ug_ceil


    if Ug>Ug_ceil:

        tmax=tmax_ceil
    else:
        tmax
        

    sol=next_state(state)
    ns.append(sol[2])
    state=sol

plt.plot(time,ns)
plt.show()


    







