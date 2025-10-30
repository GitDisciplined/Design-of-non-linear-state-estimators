import numpy as np
import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd
#import scienceplots




def hovorka(t,y,par):
    G1,G2,Q1,Q2,S1,S2,I,x1,x2,x3,C=y

    (Vg, R_thr, R_cl, V1, ka, ke, Bio, tmax,
    Ug_ceil, F01c, Fr, EGP, Ug, u, D, G,kint) = par
    
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



def meal(n):
    
    if n==100:
        
        return 500
    #elif n==400  :

        #return 100
    #elif n==700 :
        #return 50

    #elif n==800:

        #return  50

    else:

        return 0
    
    



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


Vg=20
R_thr=9
R_cl=.01
V1=.12
ka=.018
ke=.14
Bio=.8
tmax=10
temp=tmax
kint=.07
Ug_ceil=.02







#plt.style.use(['science', 'ieee'])

#plt.plot(time,ns,ns3)
#plt.show()


#plt.plot(time,ns1,ns2)
#plt.show()

#plt.plot(time,carb)
#plt.show()


#fig, axs = plt.subplots(2, 2, figsize=(6, 4))
#axs[0,0].plot(time, ns2)
#axs[0,0].set_title('Plasma glucose concentration ')
#axs[0,0].set_ylabel('Glucose concentration (mg/dL)')
#axs[0,0].set_xlabel('time (minutes)')



#axs[0,1].plot(time,ns3)
#axs[0,1].set_title('Interstitial glucose concentration')
#axs[0,1].set_ylabel('Glucose concentration (mg/dL)')
#axs[0,1].set_xlabel('time (minutes)')
#axs[1].legend()

#axs[1,0].plot(time,carb)
#axs[1,0].set_title('meal intake rate')
#axs[1,0].set_ylabel('Carb intake rate (mmol/minute)')
#axs[1,0].set_xlabel('time (minutes)')

#axs[1,1].plot(time,ns,ns1)
#axs[1,1].set_title('Glucose mass in gut and intestine')
#axs[1,1].set_ylabel('glucose mass(mmol)')
#axs[1,1].set_xlabel('time (minutes)')

#plt.tight_layout()



#plt.show()



# normal bg range 3.9-5.5
# after meal <7.8


#minimize cost
#s.t
   #constraint 1   f1(S1)
   #constarint 2   f2(G)

#target= sum over full horizon 1/no of horizons[(x2*q2)*u + (1-u)(x1*q1)]^2
   # searching u in [2,25]

          #3.9<q1/vg<7.8

    


def initialize():

    X0=np.array([0,0,90,90,0,0,0,0,0,.9,0])
    Q=np.diag([1,1,0,0,1,1,1,1,1,1,1])
    R=9
    P0=np.diag([200,100,100,100,100,100,100,100,100,100,100])

    return X0,P0,Q,R



def pos_def(a):
 

    eigen=np.linalg.eig(a)[0]
    count=0
    for ei in range(len(eigen)):
        if eigen[ei]<0:
            count=count+1
    if count==0:
        pd=0
    else:
        pd=1
    return pd



def sigma_points(dim,X,P,alpha=.001,beta=2):

    
    k=0
    la= alpha**2*(dim+k)-dim
    gamma1=math.sqrt(dim+la)
    sigma_p=[]
    sigma_n=[]

    epsilon=.00001
    temp=pos_def(P)
    while(temp==1):
        P=P+epsilon*(np.eye(dim))
        temp=pos_def(P)
        #print(P)

    
    for i in range(dim):

        
        #sigma_p.append(X+(sqrtm((dim+k)*P)[:,i]))
        #sigma_n.append(X-(sqrtm((dim+k)*P)[:,i]))

        sigma_p.append(X+(gamma1*np.linalg.cholesky(P)[:,i]))
        sigma_n.append(X-(gamma1*np.linalg.cholesky(P)[:,i]))


    #print(np.array(sigma_p))
    #print(np.array(sigma_n))
    

    return np.array(sigma_p),np.array(sigma_n),X


def weights(dim,alpha=.001,beta=2):
    k=0

    #l=alpha**2*(dim+k)-dim
    
    #weights= l/(2*(dim+l))
    #weights0_m=l/(dim+l)
    #weights0_p=(l/(dim+l)) +(1-alpha**2-beta)



    weights= 1/(2*(dim+k))
    weights0_m=k/(dim+k)
    weights0_p=(k/(dim+k))

    #print(weights)

    return weights,weights0_m,weights0_p








def predict(dim,P0,X0,Q):

 
    s_pos,s_neg,s_0=sigma_points(dim,X0,P0,alpha=.001,beta=2)
   
    com_weight,weight_mean,weight_cov=weights(dim,alpha=.001,beta=2)

    s_pos_tr=[]
    s_neg_tr=[]
    s_pos_c=[]
    s_neg_c=[]

    for m in range(len(s_pos)):
        u0_tf=next_state(s_pos[m],par)

        s_pos_tr.append(u0_tf*com_weight)
        

    for n in range(len(s_neg)):
        u_tf=next_state(s_neg[n],par)
      
        s_neg_tr.append(u_tf*com_weight)

    s_0_tr=next_state(s_0,par)*weight_mean


    pred_ns= sum( s_pos_tr)+sum(s_neg_tr)+ s_0_tr


    for o in range(len(s_pos)):

        s_pos_c.append(np.outer(s_pos_tr[o]-pred_ns,s_pos_tr[o]-pred_ns) *com_weight)

    for r in range(len(s_neg)):
        
        s_neg_c.append(np.outer(s_neg_tr[r]-pred_ns,s_neg_tr[r]-pred_ns)* com_weight)

    s_0_c=np.outer( s_0_tr-pred_ns,s_0_tr-pred_ns)*weight_cov


    pred_cov=sum( s_pos_c)+sum(s_neg_c)+ s_0_c+ Q

    #print(pred_cov)

    pred_cov=(pred_cov+pred_cov.T)/2
    
    #print(pred_ns)
    true_state=X0
    
    return pred_ns, pred_cov,true_state




def update(dim,actual_measure,pred_state,pred_covar,R):

    u_pos,u_neg, u_0=sigma_points(dim,pred_state,pred_covar,alpha=.001,beta=2)
    com_weight_u,weight_mean_u,weight_cov_u=weights(dim,alpha=.001,beta=2)

    #noise2=np.random.normal(0,R)
    H=np.array([0,0,0,0,0,0,0,0,0,0,1])
    #actual_measure=np.dot(H,actual_measure)+noise2

    u_pos_tr=[]
    u_neg_tr=[]

    for x in range(len(u_pos)):
        u_pos_tr.append(np.dot(H,u_pos[x])*com_weight_u)
    for y in range(len(u_pos)):
        u_neg_tr.append(np.dot(H,u_neg[y])*com_weight_u)


    u_0_tr= np.dot(H,u_0)*weight_mean_u


    pred_measure=sum(u_pos_tr)+sum(u_neg_tr)+u_0_tr
    #innovation_covariance

    u_pos_c=[]
    u_neg_c=[]

    for a1 in range(len(u_pos)):
        tf1=u_pos_tr[a1]
        u_pos_c.append(np.outer(tf1-pred_measure,tf1-pred_measure) *com_weight_u)
    for a2 in range(len(u_neg)):
        tf2=u_neg_tr[a2]
        u_neg_c.append(np.outer(tf2-pred_measure,tf2-pred_measure)*com_weight_u)

    u_0_c= np.outer(u_0_tr-pred_measure,u_0_tr-pred_measure)* weight_cov_u


    innov_cov=sum(u_pos_c)+sum(u_neg_c)+u_0_c+ R
   
    

   #cross_covariance

    u_pos_cc=[]
    u_neg_cc=[]

    for b1 in range(len(u_pos)):
        u_pos_cc.append((u_pos[b1]-pred_state)*(u_pos_tr[b1]-pred_measure)*com_weight_u)
    for b2 in range(len(u_neg)):
        u_neg_cc.append((u_neg[b2]-pred_state)*(u_neg_tr[b2]-pred_measure)*com_weight_u)

    u_0_cc=  (u_0-pred_state)*(u_0_tr-pred_measure)* weight_cov_u

    cross_cov=sum(u_pos_cc)+sum(u_neg_cc)+u_0_cc


    #kalman_gain

    K= cross_cov *(1/innov_cov)
   
   
   #updated_next_state
   
    updated_ns= pred_state+K*(actual_measure-pred_measure)

   #updated_covariance

    updated_cov=pred_covar-(K*innov_cov)@K.T

    updated_cov=(updated_cov+updated_cov.T)/2

    #joseph form of covariance
    #updated_cov=((K*R)@K.T)+((np.identity(3)-(K@H))@pred_covar)@(np.identity(3)-(K@H)).T

    return updated_ns, updated_cov


#def measure():
#    return np.random.normal(90,5)

#def meal_intake():

#    return np.array([300,500,0])
    



# main program

df = pd.read_csv("results.csv")
#bg= df['BG'].tolist()
ig= df['CGM'].tolist()

dim=10
X0,P0,Q,R=initialize()
state_estimate1=[]
truth_state1=[]
state_estimate2=[]
truth_state2=[]
state_estimate3=[]
truth_state3=[]
no_of_steps=len(ig)
time_steps=np.linspace(0,no_of_steps-1,no_of_steps)




for itere in range(no_of_steps):
    
    #if itere==no_of_steps/2:
      #  X0=meal_intake()

    


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


    Vg=20
    R_thr=9
    R_cl=.01
    V1=.12
    ka=.018
    ke=.14
    Bio=.8
    tmax=10
    temp=tmax
    kint=.07
    Ug_ceil=.02

    X0=X0


#initial

#X0=np.array([0,0,90,0,0,0,2,0,0,.9,0])
#u=0
#D=0
#EGP=16.9



    ns=[]
    ns1=[]
    ns2=[]
    ns3=[]
    carb=[]

    D= meal(itere)

    carb.append(D)
    

    if itere==70  :
        u=.1
    else:
        u=0
    
    Ug=X0[1]/tmax
    G= X0[2]/Vg
    

    EGP=max((1-X0[9])*EGP0,0)

    F01c=(F01s*G)/(G+1)
  
    if G>=R_thr:

        Fr=R_cl*(G-R_thr)*Vg
    else:
        Fr=0

    
    tmax_ceil=X0[1]/Ug_ceil


    if Ug>Ug_ceil:

        tmax=tmax_ceil
    else:
        tmax=temp

   

    

    par=(Vg, R_thr, R_cl, V1, ka, ke, Bio, tmax,
    Ug_ceil, F01c, Fr, EGP, Ug, u, D, G,kint)



   


    predicted_state,predicted_cov,true_state=predict(dim,P0,X0,Q)


    #print(P0)

    #print(predicted_state)

    estimated_ns, estimated_cov=update(dim,ig[itere],predicted_state,
                                       predicted_cov,R)

    #print(estimated_ns)
    #print(itere)

    

    X0,P0=estimated_ns[0], estimated_cov
    
        
    
    
    #print(X0)

    state_estimate1.append(estimated_ns[0][10])
    #truth_state1.append(true_state[2])

    state_estimate2.append(estimated_ns[0][2]/Vg)
    #truth_state2.append(true_state[1])

    #state_estimate3.append(estimated_ns[2])
    #truth_state3.append(true_state[2])
    
    #print(estimated_ns[0]-truth_state[0])
    #print(estimated_ns)


#print(state_estimate,predicted_cov)

plt.plot(time_steps,state_estimate2,state_estimate1,label='UKF')
#plt.plot(time_steps,ig,label='truth')
plt.xlabel("time")
plt.ylabel("glucose concentration")
plt.title("Estimation of blood glucose")
plt.legend()
plt.show()


#plt.plot(time_steps,state_estimate2,label='UKF')
#plt.plot(time_steps,truth_state2,label='truth')
#plt.xlabel("time")
#plt.ylabel("plasma insulin concentration")
#plt.title("Estimation of plasma insulin")
#plt.legend()
#plt.show()


#plt.plot(time_steps,state_estimate3,label='UKF')
#plt.plot(time_steps,truth_state3,label='truth')
#plt.xlabel("time")
#plt.ylabel("interstitial insulin concentration")
#plt.title("Estimation of interstitial insulin")
#plt.legend()
#plt.show()










