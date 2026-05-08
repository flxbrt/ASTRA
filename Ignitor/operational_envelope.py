# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 09:31:14 2026

@author: felix
"""



import numpy as np
import cantera as ct
from CoolProp.CoolProp import PropsSI as psi
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
import matplotlib.cm as cm


from pressure_loss_cascade import Cascade, Fixed, Resistor, Pipe



def plot_in_p_rof_space(
    p_grid, rof_grid,
    pfu_arr, pox_arr,
    dp_ref_fu_arr, dp_ref_ox_arr,
    v_fu_arr, v_ox_arr,
    Vdot_fu_arr, Vdot_ox_arr,
    fu_acc_rel_arr, ox_acc_rel_arr,
    of_min_rel_arr, of_max_rel_arr,
    dp_iso=0.20,
    Vdot_fu_limit=1000,
    Vdot_ox_limit=2000,
    p_fu_iso=22,
    p_ox_iso=22,
):
    """
    Plot all design-space quantities on a pressure × O/F grid.

    Parameters
    ----------
    p_grid, rof_grid        : 2-D meshgrids [Pa] and [-]
    pfu_arr, pox_arr        : upstream injector pressures [Pa]
    dp_ref_fu_arr, ...      : relative pressure drops [-]  (e.g. 0.20 = 20 %)
    v_fu_arr, v_ox_arr      : injection velocities [m/s]
    Vdot_fu_arr, Vdot_ox_arr: volume flow rates [std l/min]
    fu_acc_rel_arr, ...     : relative measurement uncertainties [-]
    of_min_rel_arr, ...     : relative O/F uncertainty bounds [-]
    dp_iso                  : isocontour level for relative pressure drop (default 0.20)
    Vdot_fu_limit           : isocontour for fuel volume flow [std l/min] (default 1000)
    Vdot_ox_limit           : isocontour for ox  volume flow [std l/min] (default 2000)
    """

    # ── convenience axes labels ───────────────────────────────────────────────
    p_axis   = p_grid[0, :]   / 1e5      # bar  (x-axis, one row is enough)
    rof_axis = rof_grid[:, 0]            # O/F  (y-axis, one column)

    def _labels(ax, xlabel="p_cc [bar]", ylabel="O/F [-]"):
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.tick_params(labelsize=8)

    def _pcolor(ax, data, title, unit="", cmap="viridis", fmt=None):
        """Filled contour + colorbar."""
        pc = ax.pcolormesh(p_axis, rof_axis, data, cmap=cmap, shading="auto")
        cb = plt.colorbar(pc, ax=ax, pad=0.02)
        cb.set_label(unit, fontsize=8)
        cb.ax.tick_params(labelsize=7)
        if fmt:
            cb.ax.yaxis.set_major_formatter(ticker.FuncFormatter(fmt))
        ax.set_title(title, fontsize=9, pad=4)
        _labels(ax)
        return pc

    # def _pcolor(ax, data, title, unit="", cmap="viridis", fmt=None, n_levels=None):
    #     """Filled contour + discrete colorbar."""
    #     if n_levels is None:
    #         # default: match the coarser of the two grid axes
    #         n_levels = max(p_grid.shape[1], p_grid.shape[0])
    
    #     vmin, vmax = data.min(), data.max()
    #     boundaries = np.linspace(vmin, vmax, n_levels + 1)
    #     norm = BoundaryNorm(boundaries, ncolors=256)
    
    #     pc = ax.pcolormesh(p_axis, rof_axis, data,
    #                        cmap=cmap, norm=norm, shading="auto")
    #     cb = plt.colorbar(pc, ax=ax, pad=0.02,
    #                       ticks=boundaries,
    #                       boundaries=boundaries)
    #     cb.set_label(unit, fontsize=8)
    #     cb.ax.tick_params(labelsize=7)
    #     if fmt:
    #         cb.ax.yaxis.set_major_formatter(ticker.FuncFormatter(fmt))
    #     ax.set_title(title, fontsize=9, pad=4)
    #     _labels(ax)
    #     return pc
    
    def _iso(ax, data, level, color="red", lw=1.4, label=None):
        """Single isocontour line."""
        cs = ax.contour(p_axis, rof_axis, data, levels=[level],
                        colors=[color], linewidths=lw)
        if label:
            ax.clabel(cs, fmt=lambda _: label, fontsize=7, inline=True)


    # ── figure layout  ────────────────────────────────────────────────────────
    # Rows: 1=pressure, 2=dp_ref, 3=velocity(+ratio), 4=volume flow, 5=uncertainty
    # Cols: normally 2; row 3 and 5 have a third column
    fig = plt.figure(figsize=(13, 20))
    fig.subplots_adjust(hspace=0.42, wspace=0.38)

    # GridSpec with 3 columns; rows 1-2 and 4 span 2 cols each using colspan
    from matplotlib.gridspec import GridSpec
    gs = GridSpec(5, 3, figure=fig,
                  hspace=0.42, wspace=0.38)

    # ── Row 1 : upstream pressures ────────────────────────────────────────────
    ax_pfu = fig.add_subplot(gs[0, 0])
    ax_pox = fig.add_subplot(gs[0, 1])

    # keep row 1 symmetric: both get their own half
    # ax_pfu = fig.add_subplot(gs[0, 0])
    # ax_pox = fig.add_subplot(gs[0, 1])
    # (leave gs[0,2] empty so rows 1-2 align with row 3's 3-col layout)

    _pcolor(ax_pfu, pfu_arr / 1e5,  "Max pressure in line — fuel",  "p_fu [bar]",  cmap="Blues") # cmap="plasma")
    _iso(ax_pfu, pfu_arr / 1e5, p_fu_iso, label=f"{p_fu_iso:.0f} bar")
    
    _pcolor(ax_pox, pox_arr / 1e5,  "Max pressure in line — ox",    "p_ox [bar]",  cmap="Blues") # cmap="plasma")
    _iso(ax_pox, pox_arr / 1e5, p_ox_iso, label=f"{p_ox_iso:.0f} bar")

    # ── Row 2 : relative pressure drops ──────────────────────────────────────
    ax_dpfu = fig.add_subplot(gs[1, 0])
    ax_dpox = fig.add_subplot(gs[1, 1])

    pct_fmt = lambda v, _: f"{v*100:.0f} %"

    _pcolor(ax_dpfu, dp_ref_fu_arr, "Rel. pressure drop — fuel", "Δp/p_cc [-]",
            # cmap="RdYlGn_r", fmt=pct_fmt)
            cmap="Blues", fmt=pct_fmt)
    _iso(ax_dpfu, dp_ref_fu_arr, dp_iso,
         label=f"{dp_iso*100:.0f} %")

    _pcolor(ax_dpox, dp_ref_ox_arr, "Rel. pressure drop — ox",  "Δp/p_cc [-]",
            # cmap="RdYlGn_r", fmt=pct_fmt)
            cmap="Blues", fmt=pct_fmt)
    _iso(ax_dpox, dp_ref_ox_arr, dp_iso,
         label=f"{dp_iso*100:.0f} %")

    # ── Row 3 : injection velocities + ratio ──────────────────────────────────
    ax_vfu  = fig.add_subplot(gs[2, 0])
    ax_vox  = fig.add_subplot(gs[2, 1])
    ax_vrat = fig.add_subplot(gs[2, 2])

    _pcolor(ax_vfu,  v_fu_arr,          "Injection velocity — fuel",  "v_fu [m/s]",  cmap="Blues")# cmap="coolwarm")
    _pcolor(ax_vox,  v_ox_arr,          "Injection velocity — ox",    "v_ox [m/s]",  cmap="Blues")# cmap="coolwarm")
    _pcolor(ax_vrat, v_fu_arr / v_ox_arr, "Velocity ratio v_fu / v_ox", "[-]",       cmap="Blues")# cmap="PuOr")

    # ── Row 4 : volume flow rates ─────────────────────────────────────────────
    ax_qfu = fig.add_subplot(gs[3, 0])
    ax_qox = fig.add_subplot(gs[3, 1])

    _pcolor(ax_qfu, Vdot_fu_arr, "Vol. flow — fuel",  "Vdot_fu [std l/min]", cmap="Blues")
    _iso(ax_qfu, Vdot_fu_arr, Vdot_fu_limit,
         label=f"{Vdot_fu_limit:.0f} l/min")

    _pcolor(ax_qox, Vdot_ox_arr, "Vol. flow — ox",    "Vdot_ox [std l/min]", cmap="Blues") #cmap="Oranges")
    _iso(ax_qox, Vdot_ox_arr, Vdot_ox_limit,
         label=f"{Vdot_ox_limit:.0f} l/min")

    # ── Row 5 : measurement uncertainties ────────────────────────────────────
    ax_ufu  = fig.add_subplot(gs[4, 0])
    ax_uox  = fig.add_subplot(gs[4, 1])
    ax_uof  = fig.add_subplot(gs[4, 2])

    of_worst = np.maximum(np.abs(of_min_rel_arr), np.abs(of_max_rel_arr))

    _pcolor(ax_ufu, fu_acc_rel_arr, "Fuel flow uncertainty",    "σ_fu / Vdot_fu [-]",
            # cmap="YlOrRd", fmt=pct_fmt)
            cmap="Blues", fmt=pct_fmt)
    _pcolor(ax_uox, ox_acc_rel_arr, "Ox flow uncertainty",      "σ_ox / Vdot_ox [-]",
            # cmap="YlOrRd", fmt=pct_fmt)
            cmap="Blues", fmt=pct_fmt)
    _pcolor(ax_uof, of_worst,       "Max O/F uncertainty",      "|σ_OF/OF|_max [-]",
            # cmap="YlOrRd", fmt=pct_fmt)
            cmap="Blues", fmt=pct_fmt)

    fig.suptitle("Design space — ignitor operating points", fontsize=12, y=0.995)
    plt.show()
    return fig


def ignitor_combustion(fu, ox, p, of, T_fu=293, T_ox=293):
    
    # Air: 23.2% O2, 76.8% N2 by mass (= ~21/79 by mole)
    air_o2_massfrac = 0.21
    air_n2_massfrac = 0.79
    ox_dict = {
        'O2': of * air_o2_massfrac,
        'N2': of * air_n2_massfrac
    }
        
    composition = {fu: 1.0, **ox_dict}
    
    comb = ct.Solution('h2o2.yaml')              # obtain enthalpy at equilibrium conditions
    
    comb.Y = composition
    
    T_mean = (T_fu + T_ox*of)/(1+of)        # mass averaged temperature
    
    comb.TP = T_mean, p
    comb.equilibrate('HP') 
    
    # T_ign = comb.T
    # M_ign = comb.mean_molecular_weight 
    # gamma = comb.cp_mass / comb.cv_mass  

    # c_star_id = np.sqrt(8314/M_ign*T_ign/gamma*((gamma+1)/2)**((gamma+1)/(gamma-1)))

    return comb # T_ign, M_ign, gamma, c_star_id

def ignitor_massflows(of, dth, eta, comb):
    
    Ath = np.pi*dth**2/4
    
    T_ign = comb.T
    M_ign = comb.mean_molecular_weight 
    gamma = comb.cp_mass / comb.cv_mass  

    print(T_ign, M_ign, gamma)

    c_star_id = np.sqrt(8314/M_ign*T_ign/gamma*((gamma+1)/2)**((gamma+1)/(gamma-1)))
    c_star_eff = c_star_id * eta
    
    m_tot = comb.P*Ath/c_star_eff
    
    mdot_fu = m_tot/(1+of)
    mdot_ox = m_tot*of/(1+of)
    
    return mdot_fu, mdot_ox

def convert_m_to_qn(mdot_fu, mdot_ox):
    
    p_std = 101325             # https://en.wikipedia.org/wiki/International_Standard_Atmosphere
    T_std = 288.15  
    
    rho_n_fu =  psi('D', 'P', p_std ,'T', T_std, 'H2')
    rho_n_ox =  psi('D', 'P', p_std ,'T', T_std, 'Air')
    
    return mdot_fu/rho_n_fu, mdot_ox/rho_n_ox

def injector_pressure_drop(dfu, dox, mdot_fu, mdot_ox, p, cd, T_fu=293, T_ox=293):
    
    def converge_on_pressure(mdot, A, cd, pcc, T, fluid):
        
        # pcc is static pressure
        # p is the total pressure I want to figure out
        # T is the total temperature upstream
        
        p = pcc
        T_st = T
        err = 1e5
        
        mdot = mdot/cd
        
        while err > 1:
            
            rho = psi('D', 'P', pcc ,'T', T_st, fluid)
            v = mdot / A / rho
            a = psi('speed_of_sound', 'P', pcc ,'T', T_st, fluid)
            Ma = v/a
            
            gamma = psi('Cpmass', 'P', pcc ,'T', T_st, fluid) / psi('Cvmass', 'P', pcc ,'T', T_st, fluid)
            
            p_n = pcc*(1+((gamma-1)/2)*Ma**2)**(gamma/(gamma-1))
            
            T_st = T*(1 + ((gamma-1)/2)*Ma**2)**(-1)
            
            err = abs(p_n - p)

            p = p_n
            
            # print(err)
            
        return p, v
        
    Afu = np.pi*dfu**2/4
    Aox = np.pi*dox**2/4
    
    pfu, v_fu = converge_on_pressure(mdot_fu, Afu, cd, p, T_fu, 'H2')
    pox, v_ox = converge_on_pressure(mdot_ox, Aox, cd, p, T_ox, 'Air')
    
    return pfu, pox, v_fu, v_ox
    


def max_line_pressure(mdot_fu, mdot_ox, dfu, dox, pinj_fu, pinj_ox, T_fu, T_ox, margin=0.2):
    
    fu_line = Cascade()
    fu_line.set_fluid('H2')
        
    pipe = Pipe(l=2, d=4e-3)
    filt = Resistor(kv=0.57*0.865, T=T_ox)
    mfc = Fixed(0.17e5) # 1000 slpm https://documents.alicat.com/specifications/DOC-SPECS-MCQ-HIGH.pdf
    mv = Resistor(kv=0.6, T=T_fu)
    cv = Resistor(kv=0.47*0.865, T=T_fu)
    # filt = Resistor(kv=0.5, T=T_fu)
    # inj = Fixed(dpinj_fu)
    
    fu_line.set_layout([pipe, filt, mfc, mv, cv])
    fu_line.set_boundary_condition(value=(pinj_fu, mdot_fu), typ='pm', index=(4,), port=('outlet', ))
    fu_line.solve()
    # fu_line.print_cascade()
    
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
    # ox_line.print_cascade()
    
    max_pfu = fu_line.layout[0].pin
    max_pox = ox_line.layout[0].pin
    
    return max_pfu, max_pox, max_pfu*(1+margin), max_pox*(1+margin)



def measurement_accuracy(Vdot_fu, Vdot_ox, fu_acc_abs, ox_acc_abs):
    
    p_std = 101325             # https://en.wikipedia.org/wiki/International_Standard_Atmosphere
    T_std = 288.15  
    
    rho_n_fu =  psi('D', 'P', p_std ,'T', T_std, 'H2')
    rho_n_ox =  psi('D', 'P', p_std ,'T', T_std, 'Air')
    
    of = Vdot_ox*rho_n_ox/(Vdot_fu*rho_n_fu)
    
    # conservatove estimation
    of_max = (Vdot_ox+ox_acc_abs)*rho_n_ox/((Vdot_fu-fu_acc_abs)*rho_n_fu)
    of_min = (Vdot_ox-ox_acc_abs)*rho_n_ox/((Vdot_fu+fu_acc_abs)*rho_n_fu)
    
    of_max_rel = (of_max - of)/of
    of_min_rel = (of_min - of)/of
    
    mdot_ox = Vdot_ox*rho_n_ox
    mdot_fu = Vdot_fu*rho_n_fu
    
    mdot_fu_acc_abs = fu_acc_abs*rho_n_fu
    mdot_ox_acc_abs = ox_acc_abs*rho_n_ox
    
    # Gauß error propagation https://de.wikipedia.org/wiki/Fehlerfortpflanzung
    # Assumption: errors on fuel and ox mass flow are independent
    delta_of = np.sqrt((1/mdot_fu*mdot_ox_acc_abs)**2 + (mdot_ox/mdot_fu**2*mdot_fu_acc_abs)**2)
    
    of_min = of - delta_of
    of_max = of + delta_of
    
    of_max_rel = (of_max - of)/of
    of_min_rel = (of_min - of)/of
    
    return fu_acc_abs/Vdot_fu, ox_acc_abs/Vdot_ox, of_min, of_max, of_min_rel, of_max_rel
    


def main(fu, ox, p, of, dth, eta, cd, fu_acc_abs, ox_acc_abs, T_fu, T_ox):
    
    comb = ignitor_combustion(fu, ox, p, of)
    
    mdot_fu, mdot_ox = ignitor_massflows(of, dth, eta, comb)
    
    Vdot_fu, Vdot_ox = convert_m_to_qn(mdot_fu, mdot_ox)
    
    Vdot_fu *= 60000 # conversion from std m3/s in std l/min
    Vdot_ox *= 60000
    
    pfu, pox, v_fu, v_ox = injector_pressure_drop(dfu, dox, mdot_fu, mdot_ox, p, cd)
    
    max_pfu, max_pox, max_pfu_margin, max_pox_margin = max_line_pressure(mdot_fu, mdot_ox, dfu, dox, pfu, pox, T_fu, T_ox, margin=0.2)    

    fu_acc_rel, ox_acc_rel, of_min, of_max, of_min_rel, of_max_rel = measurement_accuracy(Vdot_fu, Vdot_ox, fu_acc_abs, ox_acc_abs)
    
    return mdot_fu, mdot_ox, Vdot_fu, Vdot_ox, pfu, pox, v_fu, v_ox, fu_acc_rel, ox_acc_rel, of_min, of_max, of_min_rel, of_max_rel, max_pfu, max_pox, max_pfu_margin, max_pox_margin





#%%

if __name__ == '__main__':
    
    
    fu = 'H2'
    ox = 'Air'
    # pc = 12e5
    # dth = 0.006
    # dfu = 0.00176
    # dox = 0.005
    
    # pc = 10e5
    dth = 0.0065
    dfu = 0.00191
    dox = 0.0055
    
    eta = 0.75
    cd = 0.4
    
    T_fu = 293
    T_ox = 293
    
    ''' berücksichtigt aktuell noch nicht dass normvolumenströme nicht so übertragbar sind '''
    fu_acc_abs = 0.021*945 # in std l/min 
    ox_acc_abs = 0.021*2000 # datasheet says +/- 2% + 0.1% for analog read out

    p_range = np.linspace(2e5, 12e5, 21)
    rof_range = np.linspace(20, 81, 61)
    
    # p_range = np.array([[12e5]])
    # rof_range = np.array([[30]])
    
    p_grid, rof_grid = np.meshgrid(p_range, rof_range)
    
    # p = 12e5
    # of = 30
    
    # mdot_fu, mdot_ox, Vdot_fu, Vdot_ox, pfu, pox, v_fu, v_ox, fu_acc_rel, ox_acc_rel, of_min, of_max, of_min_rel, of_max_rel = main(fu, ox, p, of, dth, eta, cd, fu_acc_abs, ox_acc_abs, T_fu, T_ox)

    # dp_ref_fu = (pfu - p)/p
    # dp_ref_ox = (pox - p)/p
    
#%%
    shape = p_grid.shape

    
    # Pre-allocate result arrays (same shape as grids)
    mdot_fu_arr  = np.zeros(shape)
    mdot_ox_arr  = np.zeros(shape)
    Vdot_fu_arr  = np.zeros(shape)
    Vdot_ox_arr  = np.zeros(shape)
    pfu_arr      = np.zeros(shape)
    pox_arr      = np.zeros(shape)
    pinj_fu_arr      = np.zeros(shape)
    pinj_ox_arr      = np.zeros(shape)
    v_fu_arr     = np.zeros(shape)
    v_ox_arr     = np.zeros(shape)
    fu_acc_rel_arr  = np.zeros(shape)
    ox_acc_rel_arr  = np.zeros(shape)
    of_min_arr   = np.zeros(shape)
    of_max_arr   = np.zeros(shape)
    of_min_rel_arr  = np.zeros(shape)
    of_max_rel_arr  = np.zeros(shape)
    dp_ref_fu_arr   = np.zeros(shape)
    dp_ref_ox_arr   = np.zeros(shape)

    n_total = p_grid.size
    n_done = 0

    for i in range(shape[0]):
        for j in range(shape[1]):

            p   = p_grid[i, j]
            of  = rof_grid[i, j]

            n_done += 1
            print(f"[{n_done:>4d}/{n_total}]  p = {p/1e5:5.1f} bar  |  O/F = {of:.1f}  ...", end=' ')

            (mdot_fu, mdot_ox,
             Vdot_fu, Vdot_ox,
             pfu, pox,
             v_fu, v_ox,
             fu_acc_rel, ox_acc_rel,
             of_min, of_max,
             of_min_rel, of_max_rel, 
             max_pfu, max_pox, 
             max_pfu_margin, max_pox_margin) = main(fu, ox, p, of, dth, eta, cd, fu_acc_abs, ox_acc_abs, T_fu, T_ox)

            mdot_fu_arr[i, j]   = mdot_fu
            mdot_ox_arr[i, j]   = mdot_ox
            Vdot_fu_arr[i, j]   = Vdot_fu
            Vdot_ox_arr[i, j]   = Vdot_ox
            pfu_arr[i, j]       = max_pfu_margin # pfu
            pox_arr[i, j]       = max_pox_margin # pox
            pinj_fu_arr[i, j]   = pfu
            pinj_ox_arr [i, j]  = pox
            v_fu_arr[i, j]      = v_fu
            v_ox_arr[i, j]      = v_ox
            fu_acc_rel_arr[i, j]   = fu_acc_rel
            ox_acc_rel_arr[i, j]   = ox_acc_rel
            of_min_arr[i, j]    = of_min
            of_max_arr[i, j]    = of_max
            of_min_rel_arr[i, j]   = of_min_rel
            of_max_rel_arr[i, j]   = of_max_rel
            dp_ref_fu_arr[i, j] = (pfu - p) / p
            dp_ref_ox_arr[i, j] = (pox - p) / p
            
            # print(dp_ref_fu_arr, dp_ref_ox_arr)

            print(f"done  →  mdot_fu={mdot_fu*1000:.3f} g/s  mdot_ox={mdot_ox*1000:.3f} g/s  "
                  f"dp_fu={dp_ref_fu_arr[i,j]*100:.1f}%  dp_ox={dp_ref_ox_arr[i,j]*100:.1f}%")

    print(f"\nAll {n_total} design points complete.")
    
    fig = plot_in_p_rof_space(
        p_grid, rof_grid,
        pfu_arr, pox_arr,
        dp_ref_fu_arr, dp_ref_ox_arr,
        v_fu_arr, v_ox_arr,
        Vdot_fu_arr, Vdot_ox_arr,
        fu_acc_rel_arr, ox_acc_rel_arr,
        of_min_rel_arr, of_max_rel_arr,
        dp_iso=1,
        Vdot_fu_limit=945,
        Vdot_ox_limit=2000,
        p_fu_iso=22,
        p_ox_iso=22,
    )