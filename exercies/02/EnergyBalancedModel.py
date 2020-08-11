from model import Model

class EnergyBalancedModel(Model):
    def __init__(self):
        self.c = 9.96 * (10**6)
        self.e = 0.62
        self.alpha = 5.67 * (10**-8)
        self.s = 136

    def c1(self):
        return 1 / (4 * self.c)

    def c2(self):
        return (self.alpha * self.e) / self.c

    # ẏ(t) = c1*S(1 − α) − c2*y(t)^4 =: f (y(t), t)
    # Energy balance model
    def model(self,yk,tk):
        return (self.c1() * (self.s *(1 - self.alpha))) - (self.c2() * (yk**4))
