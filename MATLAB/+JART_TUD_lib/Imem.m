function I_mem = Imem(V_m, Ndisc, rvar, lvar)
    import JART_TUD_lib.const;
    % Extract and assign parameters
    d_r = (rvar - const.rdet) / (const.delta_r * const.rdet);
    d_l = (lvar - const.ldet) / (const.delta_l * const.ldet);

    Vpos = V_m >= 0;
    Vneg = V_m < 0;

    % Calculate p_X,Y parameters (ie. variability dependence)
    p1_0 = (const.p1_0_n + const.Dp1_0_r_n*d_r + const.Dp1_0_l_n*d_l) .* Vneg;
    p1_1 = (const.p1_1_n + const.Dp1_1_r_n*d_r + const.Dp1_1_l_n*d_l) .* Vneg;
    p1_2 = (const.p1_2_n + const.Dp1_2_r_n*d_r + const.Dp1_2_l_n*d_l) .* Vneg;
    p1_3 = (const.p1_3_n + const.Dp1_3_r_n*d_r + const.Dp1_3_l_n*d_l) .* Vneg;
    p1_4 = (const.p1_4_n + const.Dp1_4_r_n*d_r + const.Dp1_4_l_n*d_l) .* Vneg;

    p2_0 = const.p2_0_n * Vneg;

    p3_0 = (const.p3_0_n + const.Dp3_0_r_n*d_r + const.Dp3_0_l_n*d_l) .* Vneg;
    p3_1 = (const.p3_1_n + const.Dp3_1_r_n*d_r + const.Dp3_1_l_n*d_l) .* Vneg;


    p4_0 = const.p4_0_n * Vneg;
    p4_1 = (const.p4_1_n + const.Dp4_1_r_n*d_r + const.Dp4_1_l_n*d_l) .* Vneg;
    p4_2 = const.p4_2_n * Vneg;

    p5_0 = (const.p5_1_p + const.Dp5_1_r_p * d_r + const.Dp5_1_l_p * d_l) .* Vpos;
    p5_1 = (const.p5_1_n + const.Dp5_1_r_n*d_r + const.Dp5_1_l_n*d_l) .* Vneg + (const.p5_1_p + const.Dp5_1_r_p * d_r + const.Dp5_1_l_p * d_l) .* Vpos;
    p5_2 = (const.p5_2_n + const.Dp5_2_r_n*d_r + const.Dp5_2_r2_n*(d_r.^2)) .* Vneg + const.p5_2_p .* Vpos;
    
    p6_0 = (const.p6_0_p + const.Dp6_0_r_p * d_r) .* Vpos;
    p6_1 = const.p6_1_p * Vpos;

    p7_0 = (const.p7_0_n + const.Dp7_0_r_n*d_r + const.Dp7_0_l_n*d_l) .* Vneg + (const.p7_0_p + const.Dp7_0_r_p*d_r + const.Dp7_0_r2_p*(d_r.^2) + const.Dp7_0_l_p*d_l) .* Vpos;

    p7_1 = const.p7_1_p * Vpos;
    Dp7_2_l = (const.Dp7_2_l_p + const.Dp7_2_l_r_p*d_r + const.Dp7_2_l_r2_p*(d_r.^2)) .* Vpos;
    p7_2 = (const.p7_2_p + const.Dp7_2_r_p*d_r + const.Dp7_2_r2_p*(d_r.^2) + Dp7_2_l.*d_l) .* Vpos;
    p7_3 = const.p7_3_p * Vpos;

    p8_0 = (const.p8_0_p + const.Dp8_0_l_p*d_l) .* Vpos;
    p8_1 = const.p8_1_p * Vpos;

    p9_0 = (const.p9_0_n + const.Dp9_0_r_n*d_r + const.Dp9_0_l_n*d_l) .* Vneg;
    p9_1 = const.p9_1_n * Vneg;
    p9_2 = const.p9_2_n * Vneg;
    p9_3 = (const.p9_3_n + const.Dp9_3_r_n*d_r) .* Vneg;

    p10_0 = const.p10_0_n * Vneg + (const.p10_0_p + const.Dp10_0_l_p*d_l) .* Vpos;
    p10_1 = (const.p10_1_n + const.Dp10_1_r_n*d_r) .* Vneg + const.p10_1_p * Vpos;
    p10_2 = const.p10_2_n * Vneg + const.p10_2_p * Vpos;
    p10_3 = const.p10_3_n * Vneg;

    p11_0 = (const.p11_0_n + const.Dp11_0_r_n*d_r) .* Vneg + (const.p11_0_p + const.Dp11_0_l_p*d_l) .* Vpos;
    p11_1 = const.p11_1_n * Vneg + const.p11_1_p * Vpos;
    p11_2 = const.p11_2_n * Vneg + const.p11_2_p * Vpos;
    p11_3 = const.p11_3_n * Vneg;
    
    % Calculate p_X parameters (ie. voltage dependence)
    p1 = (p1_0.*(p1_1.*V_m + p1_2.*(V_m.^2))./(1 + p1_3.*V_m + p1_4.*(V_m.^2)) ) .* Vneg;

    p2 = p2_0 .* Vneg;

    p3 = (p3_0 + p3_1.*V_m) .* Vneg;

    p4 = (p4_0 - p4_1.*exp(-p4_2.*V_m)) .* Vneg + 1 * Vpos;

    p5 = (p5_0 + p5_1.*V_m + p5_2.*(V_m.^2)) .* Vneg + (p5_0 - p5_1.*exp(-p5_2.*V_m)) .* Vpos;

    p6 = 1 * Vneg + (p6_0 + p6_1.*V_m) .* Vpos;

    p7 = p7_0 .* Vneg + (p7_0 + p7_1.*V_m + p7_2.*exp(-p7_3.*V_m)) .* Vpos;

    p8 = 1 * Vneg + (p8_0 + p8_1.*V_m) .* Vpos;
    
    p9 =  sum([(p9_0 + (p9_1-p9_0)./(1+exp((V_m-p9_2)./p9_3))) .* Vneg, 0 * Vpos], 2, 'omitnan');

    p10 = (p10_0 + (p10_1-p10_0)./(1+exp((V_m-p10_2)./p10_3))) .* Vneg + (p10_0 + p10_1.*V_m + p10_2.*(V_m.^2)) .* Vpos;

    p11 = (1./(p11_0 + (p11_1-p11_0)./(1+exp((V_m-p11_2)./p11_3)))) .* Vneg + (p11_0 + p11_1.*V_m + p11_2.*(V_m.^2)) .* Vpos;
    
    I_mem = p1.*(p2.*(exp((log(Ndisc./const.Ndmin)-p3)./p4)-1)+(log(Ndisc./const.Ndmin)-p3));       %ExpLin part
    I_mem = I_mem + p5./(p6 + p7.*(p8.*exp(log(Ndisc./const.Ndmin)-p9)).^(-p10)).^(1./p11);   %Generalized logistic function part
    
end