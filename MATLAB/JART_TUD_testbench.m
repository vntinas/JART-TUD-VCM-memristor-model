clear;
clc;
close all;

set(groot,'defaultAxesFontSize',12)

% Load the CSV into a table
df = readtable('Sim_Results_Cadence_JART_Original.csv');

% Define simulation options
V_pos_peak = 1.3;
V_neg_peak = -1.3;
t_pos = abs(V_pos_peak);     %Note: We are considering a 1V/s slope
t_neg = abs(V_neg_peak);
t_max = 2*t_neg + 2*t_pos;
time = 0:t_max/1000:t_max;
Ndmin = 4e-3;
Ndmax = 22;
Ndinit = 10e-3;

% Create a pulse generator based on the class "pulse"
pul = JART_TUD_lib.pulse('trig',V_neg_peak,t_pos,V_pos_peak,t_neg);
V = pul.pulse_gen(time);

rd_list = [40.5e-9, 45e-9, 49.5e-9];
ld_list = [.36, .4, .44];

for i = 1:length(rd_list)
    for j = 1:length(ld_list)
        % Define target rd and ld values and the tolerance for selecting them from
        % the table
        rd = rd_list(i);
        ld = ld_list(j);
        tolerance = 1e-5;
        % Select rows based on the conditions
        selected_rows = df(abs(1 - df.ld/ld) < tolerance & abs(1 - df.rd/rd) < tolerance, :);
        
        % Define the simulation setup
        F = ode(Solver="stiff",RelativeTolerance=1e-12);
        F.ODEFcn = @(t,y) JART_TUD_lib.dNdisc_dt(pul.pulse_gen(t), y, rd, ld, Ndmin, Ndmax);
        F.InitialValue = Ndinit;
        
        % Run the simulation
        sol = solve(F,0,t_max,"Refine",100);
        
        Vsim = pul.pulse_gen(sol.Time);
        Imem = arrayfun(@(x,y) JART_TUD_lib.Imem(x,y,rd,ld), Vsim, sol.Solution);
        
        %% Plots
        fig = figure();
        fig.Color = 'white';
        tl = tiledlayout(fig,3,2);
        tl.TileSpacing = "compact";
        tl.Padding = "compact";
        title(tl,['r_d = ' num2str(rd*1e9) 'nm | l_d = ' num2str(ld) 'nm']);
        nexttile(tl,1)
        plot(time,V,'.',MarkerSize=12)
        grid on;
        grid minor;
        ax = gca;
        ax.XLim = [-.2 t_max+.2];
        xlabel('Time [s]');
        ylabel('V_M [V]');
        
        nexttile(tl,2,[3 1])
        plot(Vsim,abs(Imem),'.',MarkerSize=12);
        % plot(Vsim,abs(Imem),"-o");
        hold on;
        plot(selected_rows.Vm,abs(selected_rows.Im),LineWidth=2);
        hold off;
        ax = gca;
        ax.YScale = "log";
        ax.YLim(1) = 2e-7;
        ax.YLim(2) = 1e-3;
        grid on;
        grid minor;
        xlabel('V_M [V]');
        ylabel('I_M [A]');
        
        nexttile(tl,3)
        plot(sol.Time,Imem,'.',MarkerSize=12);
        hold on;
        plot(selected_rows.time,selected_rows.Im,LineWidth=2);
        hold off;
        grid on;
        grid minor;
        ax = gca;
        ax.XLim = [-.2 t_max+.2];
        xlabel('Time [s]');
        ylabel('I_M [A]');
        
        nexttile(tl,5)
        semilogy(sol.Time,sol.Solution,'.',MarkerSize=12);
        hold on;
        plot(selected_rows.time,selected_rows.Nd,LineWidth=2);
        hold off;
        grid on;
        grid minor;
        ax = gca;
        ax.XLim = [-.2 t_max+.2];
        xlabel('Time [s]');
        ylabel('N_d [x10^{26}m^{-3}]');
    end
end