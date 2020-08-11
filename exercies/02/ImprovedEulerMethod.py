from TimeIntegrator import TimeIntegrator

class ImprovedEulerMethod(TimeIntegrator):
    def integrate(self,tk,yk,f,delta_t):
        temp =  yk + delta_t * f(yk,tk)
        intermediate = yk + (delta_t/2) * temp
        return yk + delta_t * f(intermediate,tk)
