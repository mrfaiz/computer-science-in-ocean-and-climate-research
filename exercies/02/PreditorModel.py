from model import Model

class PreditorModel(Model):
    def __init__(self,gamma_preditor,delta_preditor,mu_preditor):
        self.gamma_preditor = gamma_preditor  #death rate
        self.delta_preditor = delta_preditor  # birth rate
        self.mu_preditor = mu_preditor

    # x prey, y preditor
    def model(self,focus_k,other_k):
        return focus_k * (self.delta_preditor  * other_k - self.gamma_preditor - self.mu_preditor * focus_k)
