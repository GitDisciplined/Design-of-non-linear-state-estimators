import numpy as np
import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
#import scienceplots






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
    
    if 100<=n<115:
        
        return 20
    #elif 200<=n<220 :

        #return 18
    #elif 400<=420:
        #return 14

    
    else:

        return 0



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
#G_basal = 5
#Q1 = G_basal * Vg
#state=np.array([0,0,Q1,.5*Q1,0,0,0,.001,.001,.001,Q1*.1])



#ns=[]
#ns1=[]

#ns2=[]
#ns3=[]
#carb=[]
#ins=[]





    

    





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


#for episodes in range(500):        

  #  l_bound=3
   # u_bound=7
   # vals=[l_bound,u_bound]
   # pr=[.5,.5]

    #P_G=np.random.choice(vals,p=pr)
    #current_state=P_G

    #epsilon=epsilon_min+(epsilon_start-epsilon_min)*np.exp(-ld*episodes)
   # epsilon=.6


#initial
G_basal = 5
Q1 = G_basal * Vg
X0=np.array([0,0,Q1,.5*Q1,0,0,0,.001,.001,.001,Q1*.1])

time=list(range(700))



gamma=.5
alpha=.6
#epsilon=.6
no_of_states=10
no_of_actions=8
epsilon_min=.01
epsilon_start=.8
ld=.001



q_matrix,state_list,action_list=initilaize(no_of_states,no_of_actions)
state_vec=[0,1*Vg,2*Vg,3*Vg,4*Vg,5*Vg,6*Vg,7*Vg,8*Vg,9*Vg]
action_vec=[10,20,30,40,100,200,300,500]





for episodes in range(100):
   
    u=0
    #instant=np.random.uniform(10,800)
    #if 10+instant<=episodes<=25+instant:
        #D= meal(np.random.normal(20,10))

    #elif 100+instant<=episodes<=125+instant:
        #D= meal(np.random.normal(20,10))


    #else:
    #D=meal(episodes)


    nstate=[]
    ins=[]

    epsilon=epsilon_min+(epsilon_start-epsilon_min)*np.exp(-ld*episodes)



  
    for i in range(len(time)):

        pres=X0
        ins.append(u)
        #if 95<=i<=100:
            #u=50

        #elif 220<=i<=225:
            #u=50

        #elif 585<=i<590:
            #u=50
        #else:
            #u=0
        if 50<=i<=65:
    
            D= np.random.normal(20,5)

        elif 250<=i<=265:

            D= np.random.normal(18,5)
            
        else:
            D=0
    

        #carb.append(D)
        #ins.append(u)
    
    
        Ug=X0[1]/tmax
        G= max(X0[2]/Vg,3)
        ig=X0[10]/Vg

        #print(G,ig)
        nstate.append(G)

        EGP=max((1-X0[9])*EGP0,0)

        F01c=(F01s*G)/(G+1)
  
        if G>=R_thr:

            Fr=R_cl*(G-R_thr)*Vg
        else:
            Fr=0


        par=(Vg, R_thr, R_cl, V1, ka, ke, Bio, tmax, kint,
        Ug_ceil, F01c, Fr, EGP, Ug, u, D, G,k12,EGP0)
        


        state=X0

        sol=next_state(X0,par)
        
            

        re=100*np.exp(-(sol[2]-X0[2]))+10*np.exp(-(X0[2]-sol[2]))

        

        X0=sol


        #if G<=3:

            #re=np.exp(-50)

        #elif G>=7:
            #re=np.exp(-10)

        #else:

            #re=np.exp(-(5-G))

        

    

        

    






        #agent



        prev_Q=q_matrix
        current_state=min(state_vec, key=lambda x: abs(x- pres[2]))
   
        current_action=action_criteria(q_matrix,action_list,state_vec.index(current_state))
        u=current_action
        
        #print(u)
        #u=0

        ns=min(state_vec, key=lambda x: abs(x- sol[2]))
        #print(pres[2],sol[2])

        reward=re
        
        updated_Q=q_update(alpha,gamma,q_matrix,state_vec.index(current_state),action_vec.index(current_action),state_vec.index(ns),reward)
    
        
        #current_state=ns

        #if np.linalg.norm(updated_Q-prev_Q)<=.01:
            #break

    q_matrix=updated_Q

   
        

    

    print(q_matrix)

    

new_q=q_matrix
policy=[]

for po in range(no_of_states):
    policy.append(action_vec[int(np.argmax(new_q[po]))] )   
    

plt.plot(time,nstate)
plt.show()

print(policy)       
        





















































