classdef const
   properties (Constant)
        % JART_TUD_params.m
        
        % Positive Imem Fitting Parameters
        p5_1_p = 1.3769e-03;
        p5_2_p = 8.1819e-02;
        Dp5_1_r_p = 2.3087e-04;
        Dp5_1_l_p = 1.3293e-07;
        p6_0_p = 1.9687e-01;
        p6_1_p = -2.1833e-02;
        Dp6_0_r_p = 2.6129e-02;
        p7_0_p = -9.7606e+01;
        p7_1_p = 7.8250e+00;
        p7_2_p = 9.9296e+01;
        p7_3_p = 7.1092e-02;
        Dp7_0_r_p = -7.4338e-01;
        Dp7_0_r2_p = 1.1713e-01;
        Dp7_2_r_p = 6.9547e-01;
        Dp7_2_r2_p = -1.3724e-01;
        Dp7_2_l_r_p = -9.0456e-03;
        Dp7_2_l_r2_p = -1.2221e-03;
        Dp7_0_l_p = 2.4377e+00;
        Dp7_2_l_p = -2.3728e+00;
        p8_0_p = 1.1713e-01;
        p8_1_p = 8.1370e-02;
        Dp8_0_l_p = -3.8320e-03;
        p10_0_p = 9.7733e-01;
        p10_1_p = 3.5214e-02;
        p10_2_p = 1.2856e-02;
        Dp10_0_l_p = 5.9623e-05;
        p11_0_p = 9.4207e-01;
        p11_1_p = 3.8953e-02;
        p11_2_p = 2.3436e-02;
        Dp11_0_l_p = -6.7239e-04;
        
        % Negative Imem Fitting Parameters
        p1_0_n = 1.1830e+00;
        p1_1_n = -2.7034e-03;
        p1_2_n = -4.5379e-06;
        p1_3_n = 9.9115e-01;
        p1_4_n = 4.4093e-01;
        Dp1_0_r_n = -6.2246e-02;
        Dp1_1_r_n = -1.8077e-04;
        Dp1_2_r_n = 1.7313e-04;
        Dp1_3_r_n = 5.7155e-03;
        Dp1_4_r_n = -9.7198e-04;
        Dp1_0_l_n = 1.1419e-01;
        Dp1_1_l_n = -2.0831e-04;
        Dp1_2_l_n = -8.9677e-05;
        Dp1_3_l_n = -2.3237e-02;
        Dp1_4_l_n = -1.8507e-03;
        p2_0_n = -2.5955e+03;
        p3_0_n = 6.8845e+00;
        p3_1_n = -5.8995e-01;
        Dp3_0_r_n = 1.2536e-01;
        Dp3_1_r_n = 6.5498e-02;
        Dp3_0_l_n = 2.5983e-01;
        Dp3_1_l_n = 8.5666e-02;
        p4_0_n = 2.5890e+03;
        p4_1_n = -2.9537e+00;
        p4_2_n = -5.4031e-01;
        Dp4_1_r_n = 8.2522e-02;
        Dp4_1_l_n = -7.2255e-02;
        p5_1_n = 6.4705e-04;
        p5_2_n = 5.1529e-05;
        Dp5_1_r_n = 1.5169e-05;
        Dp5_2_r_n = 6.7042e-07;
        Dp5_2_r2_n = 1.0756e-06;
        Dp5_1_l_n = 1.3260e-06;
        p7_0_n = 1.1708e-01;
        Dp7_0_r_n = 4.8662e-04;
        Dp7_0_l_n = 3.7351e-03;
        p9_0_n = 3.9052e+00;
        p9_1_n = 9.6130e+00;
        p9_2_n = -4.5637e-01;
        p9_3_n = 1.4310e+00;
        Dp9_0_r_n = -5.4723e-01;
        Dp9_3_r_n = 3.6000e-01;
        Dp9_0_l_n = 3.6802e-02;
        p10_0_n = 4.6925e-01;
        p10_1_n = 3.4731e+00;
        p10_2_n = -1.1871e+00;
        p10_3_n = 5.6947e-01;
        Dp10_1_r_n = 1.1444e-02;
        p11_0_n = 1.0667e+01;
        p11_1_n = 1.2812e-01;
        p11_2_n = 7.4414e-01;
        p11_3_n = 4.2381e-01;
        Dp11_0_r_n = 3.6290e-01;
        
        % Special fitting parameter (I plan to normalize the state variable in the next version and I will discard this parameter)
        Ndmin = 4e-3;
        
        % Physical constants do not change!
        e = 1.602e-19;          % elementary charge [C]
        kb = 1.3807e-23;        % Boltzman's constant  [VAs/K]
        Arichardson = 6.01e5;   % Richardson's constant [A/m^2K^2]
        mdiel = 9.10938e-31;    % electron rest mass [kg]
        h = 6.626e-34;          % Planck's constant [Js]
        zvo = 2;                % oxygen vacancy charge number
        eps_0 = 8.854e-12;      % vacuum permittivity [As/Vm]
        
        % Fitting Parameters for original model
        T0 = 293.0;          % ambient temperature [K]
        eps = 17;            % static hafnium oxide permittivity
        epsphib = 5.5;       % hafnium oxide permittivity related to image force barrier lowering
        phiBn0 = 0.18;       % nominal Schottky barrier height [eV]
        phin = 0.1;          % energy level difference between the Fermi level in the oxide and the oxide conduction band edge [eV]
        un = 4e-6;           % electron mobility [m^2/Vs]
        Ninit = 0.008;       % initial oxygen vacancy concentration in the disc [10^26/m^3]
        Nplug = 20;          % oxygen vacancy concentration in the plug [10^26/m^3]
        a = 0.25e-9;         % ion hopping distance [m]
        ny0 = 2e13;          % attempt frequency [Hz]
        dWa = 1.35;          % activation energy [eV]
        Rth0 = 15.72e6;      % thermal resistance of the Hafnium Oxide [K/W]
        rdet = 45e-9;        % radius of the filament area [m]
        lcell = 3;           % length of disc and plug region [m]
        ldet = 0.4;          % length of the disc region [m]
        Rtheff_scaling = 0.27; % scaling factor for gradual RESET
        RseriesTiOx = 650;   % series resistance of the TiOx layer [Ohm]
        R0 = 719.2437;       % line resistance for a current of 0 A [Ohm]
        Rthline = 90471.47;  % thermal resistance of the lines [W/K]
        alphaline = 3.92e-3; % temperature coefficient of the lines [1/K]
        
        % Variability-related parameters (Fixed values in [2])
        delta_r = 0.1;
        delta_l = 0.1;
   end
end
