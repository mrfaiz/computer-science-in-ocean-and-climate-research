from model import Model

class PreyModel(Model):
    def __init__(self,alpha_prey,beta_prey,lmda_prey):
        self.alpha_prey =alpha_prey  #birth rate
        self.beta_prey = beta_prey  # death rate
        self.lmda_prey = lmda_prey

    # x = prey , y = preditor
    def model(self,focus_k,other_k):
        return focus_k * (self.alpha_prey  -  self.beta_prey * other_k - self.lmda_prey * focus_k)
