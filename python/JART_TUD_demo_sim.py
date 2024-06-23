import JART_TUD_lib as vns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import solve_ivp
import pandas as pd
from matplotlib_inline.backend_inline import set_matplotlib_formats
set_matplotlib_formats('svg')


if __name__ == "__main__":
    df = pd.read_csv('Sim_Results_Cadence_JART_Original.csv')
    
    V_pos_peak = 1.3
    V_neg_peak = -1.3
    t_pos = abs(V_pos_peak)     #Note: We are considering a 1V/s slope
    t_neg = abs(V_neg_peak)
    t_max = 2*t_neg + 2*t_pos
    
    rel_tol = 1e-12
    
    ld_set = np.array([.36, .4, .44])
    # ld_set = np.array([.4,])
    rd_set = np.array([40.5e-9, 45e-9, 49.5e-9])
    # rd_set = np.array([45e-9])
    
    
    
    for ld in ld_set:
        for rd in rd_set:
            
            ######################################################################################################
            # Preparing the Figure and the Axes
            fig = plt.figure(figsize=(10, 5.5), constrained_layout=True)
            gs = gridspec.GridSpec(3, 2, figure=fig)
            # plt.rcParams['text.usetex'] = True        # <- Comment this line if you have problem with latex fonts
            plt.rcParams['font.size'] = 12
            
            ax1_1 = fig.add_subplot(gs[0, 0])
            ax1_1.grid(True, which='both')
            ax1_1.set_xlabel(r'$\mathrm{Time~[}s\mathrm{]}$', fontsize=16)
            ax1_1.set_ylabel(r'$V_\mathrm{M}~\mathrm{[V]}$', fontsize=16)

            ax2_1 = fig.add_subplot(gs[1, 0])
            ax2_1.grid(True, which='both')
            ax2_1.set_xlabel(r'$\mathrm{Time~[s]}$', fontsize=16)
            ax2_1.set_ylabel(r'$I_\mathrm{M}~\mathrm{[A]}$', fontsize=16)
            
            ax3_1 = fig.add_subplot(gs[2, 0])
            ax3_1.grid(True, which='both')
            ax3_1.set_xlabel(r'$\mathrm{Time~[s]}$', fontsize=16)
            ax3_1.set_ylabel(r'$N_\mathrm{d}~\mathrm{[A]}$', fontsize=16)
            
            ax1_2 = fig.add_subplot(gs[:, 1])
            ax1_2.grid(True, which='both')
            ax1_2.set_xlabel(r'$V_\mathrm{M}~\mathrm{[V]}$', fontsize=16)
            ax1_2.set_ylabel(r'$I_\mathrm{M}~\mathrm{[A]}$', fontsize=16)
            ######################################################################################################
            
            
            selected_rows = df[(abs(1 - df['ld']/ld) < 1e-5) & (abs(1 - df['rd']/rd) < 1e-5)]
            mem1 = vns.JART_TUD_memristor(Ninit = 0.010, lvar = ld, rvar = rd, Ndiscmin = 4e-3, Ndiscmax = 22)
            pul = vns.pulse(p_form='trig',V_tr_p=V_neg_peak,t_tr_p=t_pos,V_tr_n=V_pos_peak,t_tr_n=t_neg)
            
            sol = solve_ivp(mem1.dNdisc_dt, [0, t_max], [mem1.Ndisc], method='DOP853', args=(pul, ), rtol=rel_tol, max_step=t_max/1000)
            
            V_applied = pul.pulse_gen(sol.t)
            V_m = np.array(V_applied)
            Ndisc = np.array(sol.y[0])
            I_calc = vns.Imem(V_m, Ndisc, lvar = ld, rvar = rd)
            
            ######################################################################################################
            # Plots
            ax1_1.plot(sol.t,V_applied,'.-')
            ax1_1.plot(selected_rows['time'],selected_rows['Vm'])
            
            ax2_1.plot(sol.t,abs(I_calc),'.-')
            ax2_1.plot(selected_rows['time'],abs(selected_rows['Im']))
            ax2_1.set_yscale('log', base=10)
            ax2_1.set_ylim(1e-9)
            
            ax3_1.plot(sol.t,sol.y[0] ,'.-')
            ax3_1.plot(selected_rows['time'],selected_rows['Nd'])
            ax3_1.set_yscale('log', base=10)
            
            ax1_2.plot(V_applied,abs(I_calc),'.-', label=f'JART-TUD (Python)')
            ax1_2.plot(selected_rows['Vm'],abs(selected_rows['Im']), label=f'JART (Cadence)')
            ax1_2.set_yscale('log', base=10)
            ax1_2.set_ylim(1e-7)
            
            plt.suptitle(f'$r_d$ = {rd*1e9}nm | $l_d$ = {ld}nm')
            plt.show(block=False)
            ######################################################################################################
input("Press Enter to exit...")