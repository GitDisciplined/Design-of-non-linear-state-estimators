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
    S_Y=np.array([1,1,1])
    S_Y0=np.array([.1,.1,.1])

    return X,P,Q,R,S_Y,S_Y0


def sigma_points(dim,X,P,alpha=.001,beta=2):

    
    k=3-dim
    la= alpha**2*(dim+k)-dim
    gamma=math.sqrt(dim+la)
    for i in range(dim):

        sigma_p=[]
        sigma_n=[]
        
        sigma_p.append(X+(gamma*np.linalg.cholesky(P)[:,i]))
        sigma_n.append(X-(gamma*np.linalg.cholesky(P)[:,i]))

    return np.array(sigma_p).flatten(),np.array(sigma_n).flatten(),X


def weights(dim,alpha=.001,beta=2):

    l=alpha**2*(dim+k)-dim
    weights_m=[]
    weights_p=[]

    for j in range(dim):
        weights_m.append(l/(i+l))
        weights_p.append((l/(i+l))+(1-alpha**2-beta))

    return np.array(weights_m).flatten(),np.array(weights_p).flatten()



def predict(dim):

    s_pos=sigma_points(dim,X0,P0,alpha=.001,beta=2)[0]
    s_neg=sigma_points(dim,X0,P0,alpha=.001,beta=2)[1]
    s_0= sigma_points(dim,X0,P0,alpha=.001,beta=2)[2]
    
    mean_weight=weights(dim,alpha=.001,beta=2)[0]

    s_pos_tr=[]
    s_neg_tr=[]
    w_pos_tr=[]

    for m in range(len(s_pos)):

        s_pos_tr.append(next_state(s_pos[i]))

    for n in range(len(s_pos)):

        s_neg_tr.append(next_state(s_pos[i]))

    s_0_tr=next_state(s_0)



    
   


    



    pred_ns=
    
        
        
    return Y0




def update():



    return

        
    
    

    

        
        
    

    
    

    
