from TimeIntegrator import TimeIntegrator

class EulerMethod(TimeIntegrator):
    def integrate(self,tk,yk,f,delta_t):
        return yk + delta_t * f(yk,tk)
