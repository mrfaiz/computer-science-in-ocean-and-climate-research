import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.spatial import distance


class SpinUp:
    def __init__(self,alph,beta,gamma,delta,lmda,mu,N,delta_t,kappa,epsilon,use_improved_euler):
        self.alph = alph
        self.beta = beta
        self.gamma = gamma
        self.delta = delta
        self.lmda = lmda
        self.mu = mu 
        self.N = N
        self.kappa = kappa 
        self.epsilon = epsilon
        self.T = 100 
        self.maximum_step = int(2*N*N*kappa* self.T)
        self.delta_t = (self.T-0)/self.maximum_step
        # self.preditor_matrix = np.zeros((self.maximum_step, self.N))
        # self.prey_matrix= np.zeros((self.maximum_step, self.N))
        self.spaces = np.arange(0,100,1)
        self.diffusion_matrix = np.zeros((N,N))
        self.prey_at_final_step = np.zeros(N)
        self.preditor_at_final_step = np.ones(N)
        self.use_improved_euler = use_improved_euler


    def init_data(self):
        self.final_steps = 1
        self.delta_single_box_length = 1.0/N
        self.initialize_diffusion_matrix()
        self.initialize_preditor_prey()
        print(self.maximum_step)

       
    def initialize_preditor_prey(self):
        for i in range(0,self.N):
            self.prey_at_final_step[i] = math.sin(((i+1)*math.pi)/self.N)

    def initialize_diffusion_matrix(self):
        diffusion_factor = self.kappa/(self.delta_single_box_length * self.delta_single_box_length)
        print(diffusion_factor)
        for i in range(0,self.N):
            for j in range(0,N):
                if (i==j):
                    if (i == 0 or i == self.N-1):
                        self.diffusion_matrix[i][j] = (-1)*diffusion_factor
                    else:
                        self.diffusion_matrix[i][j] = (-2)*diffusion_factor
                else:
                    if(i == j+1 or j == i+1):
                       self.diffusion_matrix[i][j] = 1 * diffusion_factor
            

    def calculate_euclidean_norm(self,vector1,vector2):
        return distance.euclidean(vector1, vector2)

    def muliply_two_vector(self,vector1,vector2):
        return np.multiply(vector1,vector2)
    
    def matrix_multiplication(self,A,x):
       return np.matmul(A,x)
    
    # prey_f1 = matmul(D,xk) + xk * (alpha - (beta * yk) - (lambda * xk))
    def euler_method_prey(self):
        xk = self.prey_at_final_step
        change = (self.alph - (self.beta * self.preditor_at_final_step) - (self.lmda * xk))
        return self.matrix_multiplication(self.diffusion_matrix,xk) + (xk * change)

    # preditor_f2 = matmul(D,yk) + yk * ((delta * xk) - gamma - (mu * yk))   
    def euler_method_preditor(self):
        yk = self.preditor_at_final_step
        change = (self.delta * self.prey_at_final_step  - self.gamma  - (self.mu * yk))
        return self.matrix_multiplication(self.diffusion_matrix,yk) + (yk * change)


    def simulatioin(self):
        # print(self.prey_at_final_step)
        for i in range(1,self.maximum_step):
            # print(i)
            previous_index = i - 1

            prey_delta = self.euler_method_prey()
            preditor_delta = self.euler_method_preditor()

            if(self.use_improved_euler):
                intermediate_prey = self.prey_at_final_step + (dt/2.0) * prey_delta
                print(self.prey_at_final_step)
                print("************************")
                print(prey_delta)
                print("************************")
                print(intermediate_prey)
                intermediate_preditor = self.preditor_at_final_step + (dt/2.0) * preditor_delta
                change_prey = (self.alph - (self.beta * intermediate_preditor) - (self.lmda * intermediate_prey))
                change_preditor = (self.delta * intermediate_prey  - self.gamma  - (self.mu * intermediate_preditor))
                
                prey_delta = self.matrix_multiplication(self.diffusion_matrix,intermediate_prey) + (intermediate_prey * change_prey)
                # preditor_delta2 = self.matrix_multiplication(self.diffusion_matrix,intermediate_preditor) + (intermediate_preditor * change_preditor)

            prey_current = self.prey_at_final_step + self.delta_t * prey_delta
            preditor_current = self.preditor_at_final_step + (self.delta_t * preditor_delta)

            prey_distance = self.calculate_euclidean_norm(prey_current,self.prey_at_final_step)
            preditor_distance = self.calculate_euclidean_norm(preditor_current,self.preditor_at_final_step)

            # print("preditor_distance = {0} , prey = {1}".format(preditor_distance,prey_distance))
            if((prey_distance <= self.epsilon) and (preditor_distance <= self.epsilon)):
                break
            else:
                self.final_steps = self.final_steps + 1
                self.preditor_at_final_step = preditor_current
                self.prey_at_final_step = prey_current
        
        return self.final_steps


alpha = 2.0
beta = 3.0
gamma = 1.0
delta = 3.0
lmda = 1.0
mu = 1.0
kappa = 0.1
N = 100
t0 = 0.0
T = 10.0
dt = 0.01
epsilon = 10**(-4)
use_improved_euler = False

obj = SpinUp(alpha,beta,gamma,delta,lmda,mu,N,dt,kappa,epsilon,use_improved_euler)
obj.init_data()
total_steps = obj.simulatioin()

# fig,graph = plt.subplots()
# graph.set_xlabel("Space")
# graph.set_ylabel("Populations")

# graph.plot(obj.spaces,obj.preditor_at_final_step,label="Preditor")
# graph.plot(obj.spaces,obj.prey_at_final_step,label="Prey")

# fig.savefig("/home/faiz/SS_2020/Ocean/exercies/07/Graph_kappa_0_1_Improved_euler.png")

# Explicit Euler method

# For kappa = 0.1 , Final Steps required to converge 13391

