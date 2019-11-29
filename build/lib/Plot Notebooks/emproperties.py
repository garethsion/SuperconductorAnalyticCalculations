import numpy as np
from scipy.constants import constants as sc
from scipy.special import ellipk
from scipy.integrate import quad
from physicalproperties import PhysicalProperties

class EMProperties(PhysicalProperties):
    def __init__(self, a=None, b=None, d=None, pen_depth=None):
        super().__init__(a=a, b=b, d=d, pen_depth=pen_depth)
    
    def current_density(self,xlist,A=1):
        Kz = []
        for x in xlist:
            if abs(x) < self._a:
                Kz.append( 2*A / ( np.sqrt( (self._a**2 - x**2)*(self._b**2 - x**2) )) )
            elif abs(x) > self._b:
                Kz.append( -2*A / ( np.sqrt( (x**2 - self._a**2)*(x**2 - self._b**2) )) )

        # pad out Kz to make it the right length
        ls = int(.5*(len(xlist)-len(Kz)))
        klhs = Kz[0]*np.ones(ls)
        krhs = Kz[len(Kz)-1]*np.ones(ls)
        Kz = list(np.concatenate((np.concatenate((klhs,Kz)),krhs)))
        return Kz

    def current_density_two_slits(self,xlist,y=0,A=1):
        Kzx = []
        for x in xlist:
            #Phi = ((sc.mu_0 *A)/em.b) * ellipk(em.kp)
            Az0 = (-self.flux_per_length()/np.pi) * ( np.arctan((x-self.a)/abs(y)) - np.arctan((x+self.a)/abs(y)) )
            if abs(x)<=0:
                Kzx.append((2/sc.mu_0*self.pearl_length) * (self.flux_per_length()-Az0))
            if abs(x)>0: 
                Kzx.append((-2/sc.mu_0*self.pearl_length) * Az0)
        return Kzx
    
    def complex_field(self,xlist,A=1):
        Hy = []

        for i in range(len(xlist)):
            mod1 = xlist[i]**2 - self._a**2
            mod2 = self._b**2 - xlist[i]**2
            if xlist[i] < self._a:
                Hy.append( -self._a / ( (abs(xlist[i]) * np.sqrt(mod1 * mod2)) ) )
            if xlist[i] > self._a:
                Hy.append( self._a / ( (abs(xlist[i]) * np.sqrt(mod1 * mod2)) ) )
        return(Hy)

    def supercurrent_one_slit(self,xlist):
        return np.array( [(k**2)*self.qu for k in self.current_density_one_slit(xlist)] )
    
    def supercurrent(self,xlist):
        return np.array( [(k**2)*self.qu for k in self.current_density(xlist)] )
    
    def alpha(self,k,Phi,pearl_length):
        return (-Phi)/(np.pi*k*(1+k*pearl_length))
    
    def vector_potential_z_gauge(self,x,y):
        return quad(lambda k: self.alpha(k,self.flux_per_length(),self.pearl_length)*np.sin(k*x) * np.exp(-k*abs(y)), 0, np.inf) 
    
    def sgn(self,x):
        return np.sign(x) # Sign function

    def flux_per_length(self):
        return self.a*self.d*(self.flux_quantum**2)

    def f(self,x):
        return quad(lambda t: np.exp(-(x/self.pearl_length)*t)/((t**2)+1), 0, np.inf)
    
    def current_density_one_slit(self,x):
        return[( (-2*self.flux_per_length()*self.sgn(xr) ) / (np.pi * sc.mu_0 * self.pearl_length) ) * self.f(abs(xr))[0] for xr in x]
    