function MAP_exp_var = MAPexpectedvar(amplitude_list, response_bin, next_x, t_estim, s_estim, p_ab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   For MAP stepping rule variant
%   Uses the stimulus which minimized the variance of p there

amp_list = [amplitude_list, next_x];

t_min = 0.01;
t_max = 1.1;
s_min = 0.0001;
s_max = 0.8;



MAP_var = zeros(1,2);
for response = [false, true]
    resp_list = [response_bin, response];
    fun_MAP_E0 = @(t,s) pMAP(amp_list, resp_list, t, s, p_ab);
    MAP_norm = max(eps,integral2(fun_MAP_E0, t_min, t_max, s_min, s_max));
    % 1st order momentum:
    fun_MAP_E1 = @(t,s) t.*pMAP(amp_list, resp_list, t, s, p_ab);
    MAP_firstmoment = integral2(fun_MAP_E1, t_min, t_max, s_min, s_max)/MAP_norm;
    
    
    
    
    
    % variance
    fun_MAP_variance = @(t,s) ( (t - MAP_firstmoment).^2 ).*pMAP(amp_list, resp_list, t, s, p_ab);
    MAP_var(response+1) = integral2(fun_MAP_variance, t_min, t_max, s_min, s_max)/MAP_norm;
end

% Expected variance (to be minimized):
MAP_px = normcdf(next_x, t_estim, s_estim);
MAP_exp_var = (1 - MAP_px) * MAP_var(1) + MAP_px * MAP_var(2);

end
