# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 16:26:39 2026

@author: felix
"""



import numpy as np
from CoolProp.CoolProp import PropsSI as psi


class Component():
    def __init__(self, fluid=None, pin=None, pout=None, mdot=None):
        self.fluid = fluid
        self.pin = pin
        self.pout = pout
        self.mdot = mdot
        
    def density(self, p, T):
        self.rho = psi('D', 'P', p,'T', T, self.fluid)
        
    def viscosity(self, p, T):
        self.eta = psi('V', 'P', p,'T', T, self.fluid)
        
        
class Fixed(Component):
 
    def __init__(self, dp, **kwargs):
        super().__init__(**kwargs)  
        self.dp = abs(dp)
        
    def compute_pressure(self, direction):
        if direction=='forward':
            self.pout = self.pin - self.dp
        elif direction == 'backward':
            self.pin = self.pout + self.dp
 
 
class Resistor(Component):
 
    def __init__(self, kv, T, **kwargs):
        super().__init__(**kwargs)   
        self.kv = kv
        self.T = T
        
    def compute_pressure(self, direction):
        
        self.rho_n = psi('D', 'P', 1e5,'T', self.T, self.fluid)
        
        # print(self.mdot/self.rho_n)
        
        def dp_from_kv():
            # assuming subcritical conditions
            # https://www.schweizer-fn.de/stroemung/kvwert/kvwert.php
            
            if self.pout is not None:
                self.dp = abs((self.mdot/self.rho_n*3600/self.kv/514)**2)*self.rho_n*self.T/(self.pout/1e5)*1e5       
            
            elif self.pin is not None:
                err = 1
                pout = self.pin/1e5
                while err > 0.01:
                    step = 0.01
                    pout -= step
                    dp = abs((self.mdot/self.rho_n*3600/self.kv/514)**2)*self.rho_n*self.T/pout # in bar!!!
                    err = abs(dp - (self.pin/1e5 - pout))
                self.dp = dp*1e5
                
        dp_from_kv()
        
        if direction=='forward':
            self.pout = self.pin - self.dp
        elif direction == 'backward':
            self.pin = self.pout + self.dp
            
            
            
    
 
 
class Pipe(Component):
 
    def __init__(self, l, d, k=50e-6, T=293, **kwargs):
        super().__init__(**kwargs)  
        self.l = l
        self.d = d
        self.A = np.pi * d**2 / 4
        self.k = k                      # https://www.schweizer-fn.de/stroemung/rauhigkeit/rauhigkeit.php 
        self.T = T
        
    def compute_pressure(self, direction):
        
        def reynolds(vel, rho, eta):
            return vel*self.d*rho/eta

        def friction_factor(Re):
            if Re <= 2320:
                lam = 64/Re
                # print(1)
            elif Re > 2320 and Re <= 1e5:
                lam = 0.316/(Re**(0.25))
                # print(2)
            elif Re > 1e5:
            # else:
                lam_up = 10
                lam_low = 0
                err = 100
                while abs(err) > 1e-4:
                    lam_mid = np.mean([lam_low, lam_up])
                    err = 1/np.sqrt(lam_mid)+2*np.log10(2.51/Re/np.sqrt(lam_mid) + self.k/3.71/self.d)
                    if err>0:
                        lam_low = lam_mid
                    elif err<0:
                        lam_up = lam_mid
                    elif err == 0:
                        break
                lam = lam_mid
                # print(3)
            return lam
        
        dp = 0
        err = 1e5
        mdot = self.mdot 
        
        
        if self.pin is not None: 
            pout = self.pin
            
            while err > 1e2:
                pmid = (self.pin + pout) / 2
                rho = psi('D', 'P', pmid,'T', self.T, self.fluid)
                eta = psi('V', 'P', pmid,'T', self.T, self.fluid)
                vel = mdot/self.A/rho
                Re = reynolds(vel, rho, eta)
                lam = friction_factor(Re)
                dp_n = lam*self.l*rho*vel**2/(2*self.d)#/1e5
                err = abs(dp_n - dp)
                dp = dp_n
                pout = self.pin - dp
                
        elif self.pout is not None: 
            pin = self.pout
            
            while err > 1e2:
                pmid = (self.pout + pin) / 2
                rho = psi('D', 'P', pmid,'T', self.T, self.fluid)
                eta = psi('V', 'P', pmid,'T', self.T, self.fluid)
                vel = mdot/self.A/rho
                Re = reynolds(vel, rho, eta)
                lam = friction_factor(Re)
                dp_n = lam*self.l*rho*vel**2/(2*self.d)#/1e5
                err = abs(dp_n - dp)
                dp = dp_n
                pin = self.pout + dp
        
        self.dp = dp
            
        if direction=='forward':
            self.pout = self.pin - self.dp
        elif direction == 'backward':
            self.pin = self.pout + self.dp
            
        # print(self.dp/1e5)
        
        # print(vel)
        
        
        

        


class Cascade():
    
    def __init__(self, fluid=None):
        self.fluid = fluid          # fluid can be set at Network level or per-component
        self.layout = []
        self.n_comp = 0
        self.fp = None              # requires value at input port; if None, forward pass not required
        self.bp = None              # requires value at output port; if None, backward pass not required
    
    def set_fluid(self, fluid):
        self.fluid = fluid                         
                    
    def set_layout(self, layout):
        self.layout = layout
        self.n_comp = len(layout)
         
        for component in self.layout:
            if self.fluid is not None:
                # if component.fluid is None:
                component.fluid = self.fluid
            else:
                if component.fluid is None:
                    raise ValueError(
                        f'No fluid defined on Network and '
                        f'{type(component).__name__} (index '
                        f'{self.layout.index(component)}) has no fluid set.'
                    )
    
    def set_boundary_condition(self, value, typ, index, port=None):
        if typ == 'pm':
            # Apply massflow to every component (series network assumption)
            for component in self.layout:
                component.mdot = value[1]
 
            component = self.layout[index[0]]
 
            if port[0] == 'inlet':
                component.pin = value[0]
                if index[0] != 0:                      
                    self.layout[index[0] - 1].pout = value[0]
                    self.bp = index[0] - 1
                    
                self.fp = index[0]
 
            elif port[0] == 'outlet':
                component.pout = value[0]
                if index[0] != self.n_comp - 1:        
                    self.layout[index[0] + 1].pin = value[0]
                    self.fp = index[0] + 1
                self.bp = index[0]
 
        elif typ == 'pp':
            raise ValueError(f'{typ=} not supported yet')
 
        else:
            raise ValueError(f'{typ=} not supported')
 
    def solve(self):
        
        # forward pass
        if self.fp is not None: # if it is None, forward pass not required
            for idx in range(self.fp, self.n_comp):
                component = self.layout[idx]
                if component.pin is None:
                    component.pin = self.layout[idx - 1].pout
                component.compute_pressure('forward')
            
        # backward pass
        if self.bp is not None:
            for idx in range(self.bp, -1, -1):
                component = self.layout[idx]
                if component.pout is None:
                    component.pout = self.layout[idx + 1].pin
                component.compute_pressure('backward')
                
   
 
    def print_cascade(self):
        """
        Prints a formatted pressure/massflow cascade table, one row per component.
 
        Example output:
        ┌──────┬──────────────┬────────────────┬────────────────┬──────────────────┐
        │  ID  │  Type        │  p_in  [bar]   │  p_out [bar]   │  mdot  [kg/s]    │
        ├──────┼──────────────┼────────────────┼────────────────┼──────────────────┤
        │   0  │  Pipe        │        10.000  │         9.800  │          0.5000  │
        │   1  │  Fixed       │         9.800  │         9.600  │          0.5000  │
        └──────┴──────────────┴────────────────┴────────────────┴──────────────────┘
        """
        # ── column widths ──────────────────────────────────────────────────────
        col_id    = 6
        col_type  = 14
        col_p     = 16   # used for both pin and pout
        col_mdot  = 18
 
        def _fmt_p(value):
            """Convert Pa → bar and format, or show '---' if unknown."""
            if value is None:
                return '---'.center(col_p - 2)
            return f'{value / 1e5:>12.3f}'   # Pa → bar
 
        def _fmt_mdot(value):
            if value is None:
                return '---'.center(col_mdot - 2)
            return f'{value:>14.4f}'
 
        # ── border helpers ─────────────────────────────────────────────────────
        sep_id   = '─' * col_id
        sep_type = '─' * col_type
        sep_p    = '─' * col_p
        sep_mdot = '─' * col_mdot
 
        top    = f'┌{sep_id}┬{sep_type}┬{sep_p}┬{sep_p}┬{sep_mdot}┐'
        mid    = f'├{sep_id}┼{sep_type}┼{sep_p}┼{sep_p}┼{sep_mdot}┤'
        bot    = f'└{sep_id}┴{sep_type}┴{sep_p}┴{sep_p}┴{sep_mdot}┘'
 
        header = (
            f'│{"ID":^{col_id}}'
            f'│{"Type":^{col_type}}'
            f'│{"p_in  [bar]":^{col_p}}'
            f'│{"p_out [bar]":^{col_p}}'
            f'│{"mdot  [kg/s]":^{col_mdot}}│'
        )
 
        print(top)
        print(header)
        print(mid)
 
        for idx, comp in enumerate(self.layout):
            comp_type = type(comp).__name__
            row = (
                f'│{idx:^{col_id}}'
                f'│{comp_type:^{col_type}}'
                f'│{_fmt_p(comp.pin):^{col_p}}'
                f'│{_fmt_p(comp.pout):^{col_p}}'
                f'│{_fmt_mdot(comp.mdot):^{col_mdot}}│'
            )
            print(row)
 
        print(bot)
        
    def plot_cascade(self):
        pass
                


if __name__=='__main__':
    
    # fu_line = Cascade()
    # fu_line.set_fluid('H2')
    
    # d_fu = 5e-3
    
    # pipe1 = Pipe(l=2, d=d_fu)
    # mfc = Fixed(1e5)
    # mv = Resistor(kv=0.5, T=288)
    # cv = Resistor(kv=0.47*2/3*0.865, T=288)
    # filt = Resistor(kv=0.5, T=288)
    # inj = Fixed(12e5*0.2)
    
    # fu_line.set_layout([pipe1, mfc, mv, cv, filt, inj])
    # fu_line.set_boundary_condition(value=(12e5, 0.85e-3), typ='pm', index=(5,), port=('outlet', ))
    # fu_line.solve()
    # fu_line.print_cascade()
    
    
    
    
    
    
    # ox_line = Cascade()
    # ox_line.set_fluid('Air')
    
    # d_ox = 10e-3
    
    # pipe1 = Pipe(l=2, d=d_ox)
    # mfc = Fixed(1e5)
    # mv = Resistor(kv=0.54, T=288)
    # cv = Resistor(kv=1.2*0.865, T=288)
    # filt = Resistor(kv=0.484*0.869, T=288)
    # inj = Fixed(12e5*0.2)
    
    # ox_line.set_layout([pipe1, mfc, mv, cv, filt, inj])
    # ox_line.set_boundary_condition(value=(12e5, 42e-3), typ='pm', index=(5,), port=('outlet', ))
    # ox_line.solve()
    # ox_line.print_cascade()
    
    
    T_fu = T_ox = 293
    pinj_fu = pinj_ox = 14.4e5
    mdot_fu = 1e-3
    mdot_ox = 30e-3

    fu_line = Cascade()
    fu_line.set_fluid('H2')
        
    pipe = Pipe(l=2, d=6e-3)
    filt = Resistor(kv=0.88*0.865, T=T_fu)
    mfc = Fixed(0.97e5) # 1000 slpm https://documents.alicat.com/specifications/DOC-SPECS-MCQ-HIGH.pdf
    mv = Resistor(kv=0.28, T=T_fu)
    cv = Resistor(kv=0.47*0.865, T=T_fu)
    # filt = Resistor(kv=0.5, T=T_fu)
    # inj = Fixed(dpinj_fu)
    
    fu_line.set_layout([pipe, filt, mfc, mv, cv])
    fu_line.set_boundary_condition(value=(pinj_fu, mdot_fu), typ='pm', index=(4,), port=('outlet', ))
    fu_line.solve()
    fu_line.print_cascade()
    
    ox_line = Cascade()
    ox_line.set_fluid('Air')
    
    pipe = Pipe(l=2, d=10e-3)
    filt = Resistor(kv=0.88*0.865, T=T_ox)
    mfc = Fixed(0.59e5)
    mv = Resistor(kv=0.6, T=T_ox)
    cv = Resistor(kv=1.8*0.865, T=T_ox)
    # filt = Resistor(kv=0.484*0.869, T=T_ox)
    # inj = Fixed(12e5*0.2)
    
    ox_line.set_layout([pipe, filt, mfc, mv, cv])
    ox_line.set_boundary_condition(value=(pinj_ox, mdot_ox), typ='pm', index=(4,), port=('outlet', ))
    ox_line.solve()
    ox_line.print_cascade()
    
    max_pfu = fu_line.layout[0].pin * 1.2
    max_pox = ox_line.layout[0].pin * 1.2
    
    print(max_pfu)
    print(max_pox)
    
    
    
    

#%% testing

# cascade = Cascade()

# cascade.set_fluid(fluid='N2')

# random cascade 1
# cascade.set_layout([Fixed(0.2e5), Fixed(1.3e5), Resistor(kv=0.001, T=293), Fixed(2e5), Fixed(43e5)])

# random cascade 2
# cascade.set_layout([Fixed(0.2e5), Fixed(2e5), Fixed(43e5), Pipe(0.5, 1e-2), Resistor(kv=0.001, T=293)])

# resistor senity check
# cascade.set_layout([Resistor(kv=0.205*0.865, T=293)])

# pipe senity check
# cascade.set_layout([Pipe(0.5, 10e-3)])

# cascade.set_boundary_condition(value=(8.346e5, 0.01), typ='pm', index=(0,), port=('outlet', ))

# cascade.print_cascade()

# cascade.solve()

# cascade.print_cascade()