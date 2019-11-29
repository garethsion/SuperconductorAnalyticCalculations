import numpy as np
from scipy.constants import constants as sc
from scipy.special import ellipk
from scipy.integrate import quad
from physicalproperties import PhysicalProperties

class Inductance(PhysicalProperties):
    def __init__(self, a=None, b=None, d=None, pen_depth=None):
        super().__init__(a=a, b=b, d=d, pen_depth=pen_depth)
        
    def get_geometric_inductance(self):
        Lml = (sc.mu_0 * ellipk(self.kp)) / ellipk(self.k)
        for i in range(len(Lml)):
            if np.isnan(Lml[i]):
                Lml[i] = 0
        return Lml
    
    def get_Lml_approximate1(self):
        Lml1 = (sc.mu_0)/(2*np.pi) * (np.log(4*self._b/self._a) - (self._a/2*self._b)**2)

        for i in range(len(Lml1)):
            if np.isnan(Lml1[i]):
                Lml1[i] = 0
        return Lml1
    
    def get_Lml_approximate2(self):
        Lml2 = ((np.pi*sc.mu_0)/4) / np.log((8*self._a)/(self._b-self._a))

        for i in range(len(Lml2)):
            if np.isnan(Lml2[i]):
                Lml2[i] = 0
        return Lml2

    def g_clem(self,ec):
        l1 = (1/(2*(1-self.k)))
        l2 = ec*(1+self.k)
        return ( 1/( 2*(1-self.k) * ellipk(self.k)**2 ) ) * np.log( l1/l2 )
    
    def kinetic_inductance_clem_1(self,ec):
        return ( (sc.mu_0*self._pen_depth) / (2*self._a) ) * self.qu * self.g_clem(ec)
    
    def g_watanabe(self,w,s,t):
        k = w/(w+2*s)
        g = (1/(2*(1-k**2)*ellipk(k)**2))*( -np.log(t/(4*w)) - (w/(w+2*s))*np.log(t/(4*(w+2*s))) + ((2*(w+2*s))/(w+2*s)) *np.log(s/(w+s)) )
        return g
    
    def kinetic_inductance_watanabe(self,w,s):
        return sc.mu_0 * (self._pen_depth**2 /(w*self._d)) * self.g_watanabe(w,s,self._d) 