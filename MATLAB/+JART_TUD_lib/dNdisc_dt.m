function dNdisc_dt = dNdisc_dt(V_m, Ndisc, rvar, lvar, Ndiscmin, Ndiscmax)
    import JART_TUD_lib.const;
    % Check if the inputs are within the valid range
    % if any(Ndisc < Ndiscmin*(1-1e-8) | Ndisc > Ndiscmax*(1+1e-8))
    %     dNdisc_dt = NaN;
    %     return;
    % end

    % Calculate I_mem
    I_mem = JART_TUD_lib.Imem(V_m, Ndisc, rvar, lvar);

    % Calculate the series resistance and voltage
    A = pi * (rvar .^ 2);
    Rseries = const.RseriesTiOx + const.R0 * (1 + const.R0 * const.alphaline * (I_mem .^ 2) * const.Rthline);
    Vseries = I_mem .* Rseries;

    % Calculate the disc resistance and voltage
    Rdisc = lvar * 1e-9 ./ (Ndisc * 1e26 * const.zvo * const.e * const.un * A);
    Vdisc = I_mem .* Rdisc;

    % Update condition
    dNdisc_dt = zeros(size(Ndisc));

    OutOfBounds = (Ndisc < Ndiscmin & V_m > 0) | (Ndisc > Ndiscmax & V_m < 0) | Ndisc < Ndiscmin*(1-1e-4) | Ndisc > Ndiscmax*(1+1e-4);

    % if any((Ndisc < Ndiscmin & V_m > 0) | (Ndisc > Ndiscmax & V_m < 0))
    %     dNdisc_dt = Ndisc_update;
    %     return;
    % end

    % Calculate common parameters
    cvo = (const.Nplug + Ndisc) / 2 * 1e26;

    Vpos = V_m >= 0;
    Vneg = V_m < 0;

    E_ion = Vdisc ./ (lvar * 1e-9) .* Vneg + (V_m - Vseries) ./ (const.lcell * 1e-9) .* Vpos;
    Rtheff = (const.Rth0 * (const.rdet ./ rvar) .^ 2) .* Vneg + (const.Rth0 * const.Rtheff_scaling * (const.rdet ./ rvar) .^ 2) .* Vpos;
    Flim = (1 - (Ndisc ./ Ndiscmax) .^ 10) .* Vneg + (1 - (Ndiscmin ./ Ndisc) .^ 10) .* Vpos;

    % % Differentiate based on V_m sign
    % if V_m < 0
    %     E_ion = Vdisc ./ (lvar * 1e-9);
    %     Rtheff = const.Rth0 * (const.rdet / rvar) ^ 2;
    %     Flim = 1 - (Ndisc / Ndiscmax) .^ 10;
    % else
    %     E_ion = (V_m - Vseries) ./ (const.lcell * 1e-9);
    %     Rtheff = const.Rth0 * const.Rtheff_scaling * (const.rdet / rvar) ^ 2;
    %     Flim = 1 - (Ndiscmin / Ndisc) .^ 10;
    % end

    % Calculate gamma
    gamma = const.zvo * E_ion * const.a / (pi * const.dWa);

    % Calculate dWamin and dWamax
    dWamin = const.dWa * const.e * (sqrt(1 - gamma .^ 2) - gamma * pi / 2 + gamma .* asin(gamma));
    dWamax = const.dWa * const.e * (sqrt(1 - gamma .^ 2) + gamma * pi / 2 + gamma .* asin(gamma));

    % Calculate temperature
    T = I_mem .* (V_m - Vseries) .* Rtheff + const.T0;

    % Calculate ion current
    I_ion = const.zvo * const.e * cvo * const.a * const.ny0 * A .* ...
            (exp(-dWamin ./ (const.kb * T)) - exp(-dWamax ./ (const.kb * T))) .* Flim;

    % Update Ndisc
    Ndisc_update = -I_ion ./ (A * lvar * 1e-9 * const.e * const.zvo) / 1e26;
    dNdisc_dt(~OutOfBounds) = Ndisc_update(~OutOfBounds);
end
