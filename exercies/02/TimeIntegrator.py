
from abc import ABC, abstractmethod

class TimeIntegrator(ABC):
    def integrate(self,tk,yk,f,delta_t):
        pass

# print( issubclass(EulerMethod, TimeIntegrator)) 
# print( isinstance(EulerMethod(), TimeIntegrator)) 