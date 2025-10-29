import numpy as np
import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
#import scienceplots




def hovorka(t,y,par):
    G1,G2,Q1,Q2,S1,S2,I,x1,x2,x3,C=y

    (Vg, R_thr, R_cl, V1, ka, ke, Bio, tmax, kint,
    Ug_ceil, F01c, Fr, EGP, Ug, u, D, G) = par
    
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



#def meal(n):
    
 #   if n==200:
        
   #     return 100
   # elif n==400  :

   #     return 100
   # elif n==700 :
    #    return 50

   # elif n==800:

    #    return  50
    #else:

     #   return 0
    
    



#parameters

F01=11.1
F01s=F01/.85
EGP0=9
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
V1=.12
ka=.018
ke=.14
Bio=.8
tmax=5
temp=tmax
kint=.07
Ug_ceil=.02


def env(P_G,P_u):
    #parameters
    

    F01=11.1
    F01s=F01/.85
    EGP0=9
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
    V1=.12
    ka=.018
    ke=.14
    Bio=.8
    tmax=5
    temp=tmax
    kint=.07
    Ug_ceil=.02

#initial
    u=P_u
    D= np.random.normal(300,70)
    state=np.array([0,0,P_G,0,0,0,0,0,0,0,0])
#u=0
#D=0
#EGP=16.9

    ns=[]
    ns1=[]
    ns2=[]
    ns3=[]
    carb=[]
    #time=list(range(turns))


    

    

        #carb.append(D)
    

        #if i==100 or i==300 or i==600 or i==700:
            #u=.2
        #else:
            #u=0
    
    
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
   

    sol=next_state(state,par)

    if sol[2]/Vg<=3:

        re=-4

    elif sol[2]/Vg>=7:
        re=-6

    else:

        re=(5-(sol[2]/Vg))**2

    

    
   

    #print(sol[2]/Vg)
        #ns.append(sol[0])
        #ns1.append(sol[10])
        #ns2.append((sol[2]/Vg))
        #ns3.append(sol[1])
        #state=sol

    return sol[2]/Vg,re





def initilaize(no_of_states,no_of_actions):

    Q= np.zeros((no_of_states, no_of_actions))
    s=[]
    a=[]
    for i in range(no_of_states):
        s.append(i)
    for j in range(no_of_actions):
        a.append(j)
    
    return Q,s,a
    
def action_criteria(q_matrix,action_list,current_state):
    prob=[epsilon,1-epsilon]
    v=[0,1]
    temp=np.random.choice(v,p=prob)
    if temp==0:
        act=np.random.choice(action_list)
        action=action_vec[act]
    else:
        act=np.argmax(q_matrix[current_state])
        action=action_vec[act]

    return action
    
def q_update(alpha,gamma,q_matrix,
             current_state,current_action,next_state,reward):

    q_matrix[current_state][current_action]=q_matrix[current_state][current_action]+(alpha*(reward+gamma*np.max(q_matrix[next_state])-q_matrix[current_state][current_action]))
    new_Q=q_matrix

    return new_Q
    
    

#main
gamma=.5
alpha=.6
#epsilon=1
no_of_states=5
no_of_actions=7
epsilon_min=.01
epsilon_start=.6
ld=.001



q_matrix,state,action_list=initilaize(no_of_states,no_of_actions)
state_vec=[3,4,5,6,7]
action_vec=[0,.1,.2,.3,.4,.5,.6]


for episodes in range(500):        

    l_bound=70
    u_bound=30
    vals=[l_bound,u_bound]
    pr=[.5,.5]

    P_G=np.random.choice(vals,p=pr)
    current_state=P_G

    #epsilon=epsilon_min+(epsilon_start-epsilon_min)*np.exp(-ld*episodes)
    epsilon=.6
    for steps in range(1000):
        prev_Q=q_matrix
        current_state=min(state_vec, key=lambda x: abs(x- current_state))
        current_action=action_criteria(q_matrix,action_list,state_vec.index(current_state))
        ns,reward=env(current_state,current_action)
        ns= min(state_vec, key=lambda x: abs(x- ns))
        
        updated_Q=q_update(alpha,gamma,q_matrix,state_vec.index(current_state),action_vec.index(current_action),state_vec.index(ns),reward)
    
        
        
       
        current_state=ns

    if np.linalg.norm(updated_Q-prev_Q)<=.01:
            break

    q_matrix=updated_Q

    #print(q_matrix)
    
    

print(q_matrix) 

        
        





















































