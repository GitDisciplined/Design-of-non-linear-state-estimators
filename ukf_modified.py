import numpy as np
import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt





def bergman(t,y):

    G_basal=90
    I_basal=7.3
    p1=.03082
    p2=.02093
    n=.3
    p3=.00001602
    gamma=.003349
    h=89.5
    
    if (gamma*t*(y[0]-h)<0):
        dGdt = (-p1+y[2])*y[0]+p1*G_basal
        dIdt = -n*y[1]
        dXdt= -p2*y[2]+p3*(y[1]-I_basal)
    else:
        dGdt = (-p1+y[2])*y[0]+p1*G_basal
        dIdt = -n*y[1]+gamma*t*(y[0]-h)
        dXdt= -p2*y[2]+p3*(y[1]-I_basal)

    return [dGdt,dIdt,dXdt]

def next_state(y0):
    t_span = (0, 1) 

    solution = solve_ivp(bergman, t_span, y0)

    y_value = solution.y[:,-1]

    return y_value



def initialize():

    X=np.array([[0,0,0]])
    Q=np.diag([.1,.2,.3])
    R=np.diag([.5,.5,.5])
    P=np.diag([.1,.1,.1])

    return X,P,Q,R


def sigma_points(dim,X,P,alpha=.001,beta=2):

    
    k=3-dim
    la= alpha**2*(dim+k)-dim
    gamma=math.sqrt(dim+la)
    for i in range(dim):

        sigma_p=[]
        sigma_n=[]
        
        sigma_p.append(X+(gamma*np.linalg.cholesky(P)[:,i]))
        sigma_n.append(X-(gamma*np.linalg.cholesky(P)[:,i]))

    return np.array(sigma_p),np.array(sigma_n),X


def weights(dim,alpha=.001,beta=2):

    l=alpha**2*(dim+k)-dim
    
    weights= l/(2*(dim+l))
    weights0_m=l/(dim+l)
    weights0_p=(l/(dim+l))+(1-alpha**2-beta)

    return weights,weights0_m,weights0_p



def predict(dim):
    
    X0,Q,R,P0=initialize()
    
    s_pos=sigma_points(dim,X0,P0,alpha=.001,beta=2)[0]
    s_neg=sigma_points(dim,X0,P0,alpha=.001,beta=2)[1]
    s_0= sigma_points(dim,X0,P0,alpha=.001,beta=2)[2]
    
    mean_weight=weights(dim,alpha=.001,beta=2)[0]

    s_pos_tr=[]
    s_neg_tr=[]
    s_pos_c=[]
    s_neg_c=[]
   

    for m in range(len(s_pos)):

        s_pos_tr.append(next_state(s_pos[m])*weights(dim,alpha=.001,beta=2)[0])
        

    for n in range(len(s_neg)):

        s_neg_tr.append(next_state(s_neg[n])*weights(dim,alpha=.001,beta=2)[0])

    s_0_tr=next_state(s_0)*weights(dim,alpha=.001,beta=2)[1]


    pred_ns= sum( s_pos_tr)+sum(s_neg_tr)+ s_0_tr


    for 0 in range(len(s_pos)):

        s_pos_c.append(np.outer(next_state(s_pos_tr[o])-pred_ns,next_state(s_pos_tr[o])-pred_ns) *
                       weights(dim,alpha=.001,beta=2)[0])

    for r in range(len(s_neg)):

        s_neg_c.append(np.outer(next_state(s_neg_tr[r])-pred_ns,next_state(s_neg_tr[r])-pred_ns)*
                       weights(dim,alpha=.001,beta=2)[0])

    s_0_c=np.outer(next_state(s_0)-pred_ns)*weights(dim,alpha=.001,beta=2)[2]


    pred_cov=sum( s_pos_c)+sum(s_neg_c)+ s_0_c+ Q


    
    return pred_ns, pred_cov




def update(dim,actual_measure):

    pred_state=predict(dim)[0]
    pred_covar=predict(dim)[1]
    
    X0,Q,R,P0=initialize()

    u_pos=sigma_points(dim,X0,P0,alpha=.001,beta=2)[0]
    u_neg=sigma_points(dim,X0,P0,alpha=.001,beta=2)[1]
    u_0= sigma_points(dim,X0,P0,alpha=.001,beta=2)[2]

    H=np.array([1,0,0])

    u_pos_tr=[]
    u_neg_tr=[]

    for x in range(len(u_pos)):
        u_pos_tr.append(np.dot(H,u_pos[x])*weights(dim,alpha=.001,beta=2)[0])
    for y in range(len(u_pos)):
        u_neg_tr.append(np.dot(H,u_neg[y])*weights(dim,alpha=.001,beta=2)[0])


    u_0_tr= np.dot(H,u_0)*weights(dim,alpha=.001,beta=2)[1]


    pred_measure=sum(u_pos_tr)+sum(u_neg_tr)+u_0_tr


    

   #innovation_covariance

    u_pos_c=[]
    u_neg_c=[]

   for a1 in range(len(u_pos)):
        u_pos_c.append(np.outer(u_pos_tr[a1]-pred_measure,u_pos_tr[a1]-pred_measure)*weights(dim,alpha=.001,beta=2)[0])
   for a2 in range(len(u_neg)):
        u_neg_c.append(np.outer(u_neg_tr[a2]-pred_measure,u_neg_tr[a2]-pred_measure)*weights(dim,alpha=.001,beta=2)[0])

   u_0_c= np.outer(u_0_tr-pred_measure,u_0_tr-pred_measure)* weights(dim,alpha=.001,beta=2)[2])

   innov_cov=sum(u_pos_c)+sum(u_neg_c)+u_0_c+ R
   
    

   #cross_covariance

   u_pos_cc=[]
   u_neg_cc=[]

   for b1 in range(len(u_pos)):
        u_pos_cc.append(np.matmul(u_pos[b1]-pred_state,u_pos_tr[b1]-pred_measure)*weights(dim,alpha=.001,beta=2)[0])
   for b2 in range(len(u_neg)):
        u_neg_cc.append(np.matmul(u_neg[b1]-pred_state,u_neg_tr[b1]-pred_measure)*weights(dim,alpha=.001,beta=2)[0])

   u_0_cc=  np.matmul(u_0-pred_state,u_0_tr-pred_measure)* weights(dim,alpha=.001,beta=2)[2])

   cross_cov=sum(u_pos_cc)+sum(u_neg_cc)+u_0_cc


    #kalman_gain

   K= np.matmul(cross_cov,np.linalg.inv(innov_cov))
   
   
   #updated_next_state
   
    updated_ns= pred_state+k*(actual_measure-pred_measure)

   #updated_covariance

    updated_cov=pred_covar-K*innov_cov*innov_cov.T

    

    

   return updated_ns, updated_cov



def measure():
     return np.random.normal(2,1.4)




        
dim=3

for itere in range(100):

    

    estimated_ns, estimated_cov=update(dim,measure)

    state_estimate=[]
    
    state_estimate.append(updated_ns)


print(state_estimate)

    
    

    

        
        
    

    
    

    
