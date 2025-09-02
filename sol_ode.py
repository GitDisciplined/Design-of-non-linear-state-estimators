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
        dIdt = -n*y[1]+np.random.normal(2,5)
        dXdt= -p2*y[2]+p3*(y[1]-I_basal)+np.random.normal(1,5)
    else:
        dGdt = (-p1+y[2])*y[0]+p1*G_basal+np.random.normal(1,5)
        dIdt = -n*y[1]+gamma*t*(y[0]-h)+np.random.normal(1,5)
        dXdt= -p2*y[2]+p3*(y[1]-I_basal)+np.random.normal(1,5)

    return [dGdt,dIdt,dXdt]



y0 = [1,1,1]

t_span = (0, 10) # Solve from t=0 to t=10
t_eval = np.linspace(t_span[0], t_span[1], 100) # Evaluate solution at 100 points



    
solution = solve_ivp(bergman, t_span, y0,t_eval=t_eval)


time = solution.t
y_values = solution.y

#print(y_values)

# Plot the results
plt.plot(time, y_values[0], label='G(t) versus time')
plt.plot(time, y_values[1], label='I(t) versus time')
plt.plot(time, y_values[2], label='X(t) versus time')


plt.xlabel('Time')
plt.ylabel('Glucose/Insulin Concentration')
plt.title('Solution of ODE System')
plt.legend()
plt.grid(True)
plt.show()
