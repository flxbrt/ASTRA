 # -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
# from matplotlib import colors, colorbar, cm


class SPL_NASA():
    # NASA1971_Acoustic_Loads_generated_by_the_Propulsion_System
    
    def __init__(self):      
        print('initalized')
        
    def __call__(self, thrust, engines, velocity, diameter, f_list, distance, angle, reduction = True):
        splbp_list = list()
        
        if reduction:
            reduction_list = self.a_classification(f_list)
        else:
            reduction_list = [0]*len(f_list)
        
        for f_index, f in enumerate(f_list):
            
            oap_val = self.oap(thrust, engines, velocity)
            
            ospl_val = self.ospl(oap_val)
            
            self.strouhal(f, diameter, velocity)
            
            cabw_val = self.cabw(velocity, diameter, ospl_val, f)
            
            splbp_list.append(self.splbp(cabw_val, distance, angle, reduction_list[f_index]))
            
            osplp_val = self.osplp(splbp_list)
        
        
        
        return splbp_list, osplp_val#, ospl_val
        
    
    def a_classification(self, f_list):
        freq = [20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000]
        a_class = [-50, -30, -20, -10, -4, 0, 2, 2, -4, -9]
        
        reduction = []
        for f in f_list:
            reduction.append(np.interp(f, freq, a_class))
        return reduction
    
    def oap(self, thrust, engines, velocity, eta=1e-2): # overall acoustic power [W]
        # eta = acoustic efficiency --> 1% is conservative already
        return 0.5*eta*engines*thrust*velocity
        
    def ospl(self, oap_val): # overall sound power level [dB]
        return 10*np.log10(oap_val) + 120
        # return 10*np.log10(oap_val/1e-12) # test if log in literature means log to the base 10 or e --> it's log to the base 10
    
    def strouhal(self, f, diameter, velocity):
        self.st = f*diameter/velocity
        # return f*diameter/velocity
    
    def allocate_nrspsl(self, st_val): # normalized relative sound power spectrum level [dB]
        st = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0] # for hydrogen combustion
        # st = [0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0] # for non hydrogen combustion
        if st_val < st[0] or st_val > st[-1]:
            print('Warning, st_val out of domain of definition [0.002; 5]\n' +
                  'Return value is therefore invalid')
        # nrspsl = [0, 6.5, 9.5, 11, 8, 6, 0, -7, -13, -20, -27] # for hydrogen combustion
        nrspsl = [13, 15, 15, 13, 11.5, 8, 6, 0, -7, -13, -20, -27] # for non hydrogen combustion
        return np.interp(st_val, st, nrspsl)
    
    def cabw(self, velocity, diameter, ospl_val, f): # conventional acoustic bandwidth / sound power level in the band centered on frequency b [dB]
        # return 10*np.log10(ap/oap_val*velocity/diameter)+ospl_val-10*np.log10(velocity/diameter)+10*np.log10(delta_f)   
        delta_f = f/np.sqrt(2) # see Deutsch2019_Design_of_a_Rocket_Engine_Testbench p. 20
        return self.allocate_nrspsl(self.st)+ospl_val-10*np.log10(velocity/diameter)+10*np.log10(delta_f)   
    
    
    def get_DI(self, st, angle): # directivity index [dB]
        # round up st if angle > 60° and round down st if angle < 60° --> more conservative (see figure p. 12)
        def st_0004(angle):
            theta = [180, 160, 140, 120, 100, 80, 60, 40, 30, 20]
            theta.reverse() # reverse is required since numpy can only interpolate in case of increasing numbers
            DI = [-17.2, -15.7, -14, -11.7, -9, -4.7, 0.6, 5.9, 8.1, 4.1]
            DI.reverse()
            return np.interp(angle, theta, DI)
        def st_00125(angle):
            theta = [180, 160, 140, 120, 100, 80, 60, 40, 20]
            theta.reverse()
            DI = [-16.6, -15.1, -13.3, -11, -8, -3.6, 2.5, 6.1, 0.2]
            DI.reverse()
            return np.interp(angle, theta, DI)
        def st_004(angle):
            theta = [180, 160, 140, 120, 100, 80, 60, 47, 40, 20]
            theta.reverse()
            DI = [-15.6, -13.9, -11.9, -9.7, -6.5, -2.1, 3.3, 5.6, 5.1, -0.5]
            DI.reverse()
            return np.interp(angle, theta, DI)
        def st_0125(angle):
            theta = [180, 160, 140, 120, 100, 80, 60, 50, 40, 20]
            theta.reverse()
            DI = [-14.4, -12.7, -10.8, -8.4, -5.1, -0.7, 4, 5.1, 4.1, -1.4]
            DI.reverse()
            return np.interp(angle, theta, DI)
        def st_04(angle):
            theta = [180, 160, 140, 120, 100, 80, 60, 54, 40, 20]
            theta.reverse()
            DI = [-13.1, -11.5, -9.5, -7, -3.8, 1, 4.3, 4.5, 3.2, -2.1]
            DI.reverse()
            return np.interp(angle, theta, DI)
        
        if angle > 50:
            if st <= 0.004:
                return st_0004(angle)
            elif st <= 0.0125 and st > 0.004:
                return st_00125(angle)
            elif st <= 0.04 and st > 0.0125:
                return st_004(angle)
            elif st <= 0.125 and st > 0.04:
                return st_0125(angle)
            elif st <= 0.4 and st > 0.125:
                return st_04(angle)
            else: 
                print('Strouhal out of valid range')
                return 10 # conservative default value
        elif angle <= 50 and angle >=20:
            if st > 0.4:
                return st_04(angle)
            elif st > 0.125 and st <= 0.4:
                return st_00125(angle)
            elif st <= 0.04 and st <= 0.125:
                return st_004(angle)
            elif st <= 0.0125 and st <= 0.04:
                return st_0125(angle)
            elif st <= 0.004 and st <= 0.0125:
                return st_04(angle)
            else: 
                print('Strouhal out of valid range')
                return 10 # conservative default value
        else:
            print('Angle out of valid range')
            return 10 # conservative default value
        
    def splbp(self, cabw_val, distance, angle, red_val): # sound pressure level around band f at point p [dB]
        # print(red_val)    
        return cabw_val-10*np.log10(distance**2) - 11 + self.get_DI(self.st, angle) + red_val
    
    def osplp(self, splbp_list): # overall sound pressure level at point p [dB]
        sum_antilog = 0
        for ii in splbp_list:
            sum_antilog += 10**(ii/10)
        return 10*np.log10(sum_antilog)







class SPL_ADP(): # probably not suitable, since correlation for sound power is proportional to the exit velocity to the power of 8 --> this is commonly used for subsonic flow
    
    def __init__(self):
        print('initalized')
    
    
    def __call__(self):
        P = self.sp()
        spl_val = self.spl(P)
        print(spl_val)
        
    
    def sp(self): # sound power [W]
        rho_0 = 1 # ambient air density
        a_0 = 343.2 # ambient air speed of sound
        S = 0.1**2/4*np.pi # cross sectional area at the engine outlet
        velocity = 4500 # outlet velocity
        rho_S = 0.1 # density of the jet 
        omega = 2 # constant from table
        F = 10**0.14 # constant from table
        
        return 6.67e-5*rho_0*a_0**3*S*(velocity/a_0)**8*(rho_S/rho_0)**omega*F
    
    
    def spl(self, P): # sound power level [dB]
        return 10*np.log10(P/1e-12)
    
    



if __name__ == '__main__':

    thrust = 60
    velocity = 700
    
    dist = 5
        
    spl_nasa = SPL_NASA()
    
    oap_val = spl_nasa.oap(thrust, 1, velocity)
    
    print(spl_nasa.ospl(oap_val))