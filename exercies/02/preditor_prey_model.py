# ẋ = x(α − βy − λx)
# ẏ = y (δx − γ − μy ),


class PreditorPreyModel():
    def __init__(self, delta_t):
        self.delta_t = delta_t
        self.prey = [100]
        self.preditor = [10]
        self.alpha_prey = 0.5  #birth rate
        self.beta_prey = 0.015  # death rate
        self.gamma_preditor = 0.5  #death rate
        self.delta_preditor = 0.015  # birth rate
        self.lmda_prey = 0.0
        self.mu_preditor = 0.0

    def intial_prey_value(self, y0):
        self.prey = [y0]

    def intial_preditor_value(self, y0):
        self.preditor = [y0]

    def temp_model_function(self):
        return 1

    def time_integrator(self, k, model_function):
        preditor_model = (
            -self.gamma_preditor * self.preditor[k] +
            self.delta_preditor * self.preditor[k] * self.prey[k])
        prey_model = (self.alpha_prey * self.prey[k] -
                      self.beta_prey * self.preditor[k] * self.prey[k])
        yk_plus_1_prey = self.prey[k] + self.delta_t * prey_model
        yk_plus_1_preditor = self.preditor[k] + self.delta_t * preditor_model
        self.prey.append(yk_plus_1_prey)
        self.preditor.append(yk_plus_1_preditor)

    def get_prey(self):
        return self.prey

    def get_preditor(self):
        return self.preditor
