## Implementation of glucose and insulin estimation algorithm using unsceneted kalman filter

The pancreas is a vital organ that helps to maintain glucose level in the blood stream. When the glucose level is high, the islet cells in the pancreas produce insulin which helps to maintain the optimal level of glucose and enhances cellular absorption of glucose  for energy production.
For people with diabetes mellitus, the pancreas no longer produces insulin which can lead to both hyperglycemia and hypoglycemia. In this case, we need to continuously monitor the glucose level and determine the amount of insulin to be injected in the blood. The continuous glucose monitoring is done via a sensor that returns concentration of glucose. The automatic delivery system (ADS) uses this information to optimally control the insulin pump to inject insulin in the blood stream. The effectiveness of the decision taken by the ADS is dependant on the algorithm implemented on it. Here, algorithmic aspects of automatic insulin delivery system  is focused, especially the estimation part of the algorithms and also
 an attempt to solve Bergmans model to estimate the concentration of glucose  using Runge-kutta method and unscented kalman filter is briefly explained





![n1](https://github.com/user-attachments/assets/058433f9-0226-4264-bc7d-ce55d7eac8c2)
![n2](https://github.com/user-attachments/assets/830405d9-e1aa-4f20-87f3-1f9f13f87adf)
![n3](https://github.com/user-attachments/assets/f00a0e08-4651-4699-b2bc-41d2c5386674)




## Math for UKF

![1](https://github.com/user-attachments/assets/0eda2916-b85e-4fae-b0b5-61d441196691)

![2](https://github.com/user-attachments/assets/5d0046cc-de12-46e7-af13-4ecebe0b7994)

![3](https://github.com/user-attachments/assets/ace053f2-b49a-48a9-a3ed-cadff0c481f7)

![4](https://github.com/user-attachments/assets/dad55a04-1aa6-4b9f-93e7-09fd4a054cb9)








## Results

![glucose_estimation](https://github.com/user-attachments/assets/2c80e9ea-1452-414f-8ffd-52200a44f796)

![plasma_insulin_estimation](https://github.com/user-attachments/assets/9e218c31-43ba-402e-abb3-14961a7a44d6)

![interstitial_plasma_estimation](https://github.com/user-attachments/assets/d8b35643-2a22-430c-b6bb-96e764711e51)






