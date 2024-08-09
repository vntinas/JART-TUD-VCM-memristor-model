classdef pulse
    properties
        p_form
        V_dc
        V_tr_p
        t_tr_p
        V_tr_n
        t_tr_n
        V_sq_low
        V_sq_high
        t_sq_d
        t_sq_r
        t_sq_width
        t_sq_f
        t_sq_wait
        V_sin_A
        t_sin_per
        phi_0
        V_sin_offset
    end
    
    methods
        function obj = pulse(p_form, varargin)
            obj.p_form = p_form;
            if strcmp(p_form, 'DC')
                obj.V_dc = varargin{1};
            elseif strcmp(p_form, 'trig')
                obj.V_tr_p = varargin{1};
                obj.t_tr_p = varargin{2};
                obj.V_tr_n = varargin{3};
                obj.t_tr_n = varargin{4};
            elseif strcmp(p_form, 'square')
                obj.V_sq_low = varargin{1};
                obj.V_sq_high = varargin{2};
                obj.t_sq_d = varargin{3};
                obj.t_sq_r = varargin{4};
                obj.t_sq_width = varargin{5};
                obj.t_sq_f = varargin{6};
                obj.t_sq_wait = varargin{7};
            elseif strcmp(p_form, 'sin')
                obj.V_sin_A = varargin{1};
                obj.t_sin_per = varargin{2};
                obj.phi_0 = varargin{3};
                obj.V_sin_offset = varargin{4};
            end
        end
        
        function V = pulse_gen(obj, t)
            % Generate the pulse waveform based on the given time t
            V = zeros(size(t));
            if strcmp(obj.p_form, 'DC')
                V = obj.V_dc * ones(size(t));
            elseif strcmp(obj.p_form, 'trig')
                for i = 1:length(t)
                    t_mod = mod(t(i), 2*obj.t_tr_p + 2*obj.t_tr_n);
                    if t_mod <= obj.t_tr_p
                        V(i) = obj.V_tr_p * (t_mod / obj.t_tr_p);
                    elseif t_mod <= 2*obj.t_tr_p
                        V(i) = 2 * obj.V_tr_p - obj.V_tr_p * (t_mod / obj.t_tr_p);
                    elseif t_mod <= 2*obj.t_tr_p + obj.t_tr_n
                        V(i) = obj.V_tr_n * ((t_mod - 2*obj.t_tr_p) / obj.t_tr_n);
                    else
                        V(i) = 2 * obj.V_tr_n - obj.V_tr_n * ((t_mod - 2*obj.t_tr_p) / obj.t_tr_n);
                    end
                end
            elseif strcmp(obj.p_form, 'square')
                for i = 1:length(t)
                    t_mod = mod(t(i), obj.t_sq_d + obj.t_sq_r + obj.t_sq_width + obj.t_sq_f + obj.t_sq_wait);
                    if t_mod <= obj.t_sq_d
                        V(i) = obj.V_sq_low;
                    elseif t_mod <= obj.t_sq_d + obj.t_sq_r
                        V(i) = obj.V_sq_low + (obj.V_sq_high - obj.V_sq_low) * (t_mod - obj.t_sq_d) / obj.t_sq_r;
                    elseif t_mod <= obj.t_sq_d + obj.t_sq_r + obj.t_sq_width
                        V(i) = obj.V_sq_high;
                    elseif t_mod <= obj.t_sq_d + obj.t_sq_r + obj.t_sq_width + obj.t_sq_f
                        V(i) = obj.V_sq_high - (obj.V_sq_high - obj.V_sq_low) * (t_mod - obj.t_sq_d - obj.t_sq_r - obj.t_sq_width) / obj.t_sq_f;
                    else
                        V(i) = obj.V_sq_low;
                    end
                end
            elseif strcmp(obj.p_form, 'sin')
                for i = 1:length(t)
                    V(i) = obj.V_sin_A * sin(2 * pi * t(i) / obj.t_sin_per + obj.phi_0) + obj.V_sin_offset;
                end
            end
        end
    end
end