import numpy as np
from scipy.constants import constants as sc

class PhysicalProperties(object):
    def __init__(self,a=None,b=None,d=None,pen_depth=None):
        self._a = a
        self._b = b
        self._d = d
        self._pen_depth = pen_depth
        
        return 
    
    @property
    def a(self): 
        return self._a
    
    @a.setter 
    def a(self, a): 
        self._a = a 
        
    @property
    def b(self): 
        return self._b
    
    @b.setter 
    def b(self, b): 
        self._b = b 
        
    @property
    def d(self): 
        return self._d
    
    @d.setter 
    def d(self, d): 
        self._d = d
    
    @property
    def pen_depth(self): 
        return self._pen_depth
    
    @pen_depth.setter 
    def pen_depth(self, pen_depth): 
        self._pen_depth = pen_depth 
        
    # Modulus of complete elliptic integral of first kind
    @property
    def k(self): 
        return self._a / self._b
    
    # Complementary modulus of complete elliptic integral of first kind
    @property
    def kp(self): 
        return np.sqrt(1-self.k**2)
    
    @property
    def qu(self):
        u = self._d / self._pen_depth
        return (np.sinh(u) + u) / (8*np.sinh(u/2)**2)
    
    @property
    def pearl_length(self):
        return (2*self._pen_depth**2)/self._d

    @property
    def flux_quantum(self):
        return sc.h / (2*sc.e)

    @property
    def flux_per_length(self):
        return self._a*self._d*(self.flux_quantum()**2)
    