classdef JART_TUD_memristor
    properties
        Ndisc
        rvar
        lvar
        Ndiscmin
        Ndiscmax
    end
    
    methods
        function obj = JART_TUD_memristor(Ninit, rvar, lvar, Ndiscmin, Ndiscmax)
            obj.Ndisc = Ninit;
            obj.rvar = rvar;
            obj.lvar = lvar;
            obj.Ndiscmin = Ndiscmin;
            obj.Ndiscmax = Ndiscmax;
        end
        
        function dNdisc = dNdisc_dt(obj, t, pulse)
            dNdisc = JART_TUD_lib.dNdisc_dt(pulse.pulse_gen(t), objNdisc, obj.rvar, obj.lvar, obj.Ndiscmin, obj.Ndiscmax);            
        end
        
        function I_mem = Imem(obj, V_m)
            I_mem = JART_TUD_lib.Imem(V_m, obj.Ndisc, obj.rvar, obj.lvar);
        end
    end
end
