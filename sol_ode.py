from scipy.integrate import solve_ivp
import numpy as np
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
        dGdt = (-p1+y[2])*y[0]+p1*G_basal+np.random.normal(1,5)
        dIdt = -n*y[1]+np.random.normal(1,5)
        dXdt= -p2*y[2]+p3*(y[1]-I_basal)+np.random.normal(1,5)
    else:
        dGdt = (-p1+y[2])*y[0]+p1*G_basal+np.random.normal(1,5)
        dIdt = -n*y[1]+gamma*t*(y[0]-h)+np.random.normal(1,5)
        dXdt= -p2*y[2]+p3*(y[1]-I_basal)+np.random.normal(1,5)

    return [dGdt,dIdt,dXdt]


y0 = [1,1,1]



t_span = (0, 1) 

solution = solve_ivp(bergman, t_span, y0)

y_value = solution.y[:,-1]

print(y_value)


