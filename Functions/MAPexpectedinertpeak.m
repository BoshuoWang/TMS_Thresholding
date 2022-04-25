function MAP_exp_inert = MAPexpectedinertpeak(amplitude_list, response_bin, next_x, t_estim, s_estim, p_ab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   For MAP stepping rule variant
%   Uses the stimulus which minimized the expected momentum with respect to peak of p there 

amp_list = [amplitude_list, next_x];

t_min = 0.01;
t_max = 1.1;
s_min = 0.0001;
s_max = 0.8;

MAP_options = optimset('Display','off','MaxFunEvals',50000,'FunValCheck','on','MaxIter',50000,'TolX',1e-10);

MAP_inert = zeros(1,2);
for response = [false, true]
    resp_list = [response_bin, response];
    fun_MAP_E0 = @(t,s) pMAP(amp_list, resp_list, t, s, p_ab);
    MAP_norm = max(eps,integral2(fun_MAP_E0, t_min, t_max, s_min, s_max));
    % MAP estimation
    fun_MAP_min = @(theta) -1*logMAP(amp_list, resp_list, theta(1), theta(2), p_ab);
    [MAP_estimation, ~, map_exitflag] = fminsearch(fun_MAP_min, [t_estim, s_estim], MAP_options);
    if map_exitflag
        MAP_t_estim = MAP_estimation(1);
    else        % use current:
        MAP_t_estim = t_estim;
    end
    % inertia
    fun_MAP_InertPeak = @(t,s) ( (t - MAP_t_estim).^2).*pMAP(amp_list, resp_list, t, s, p_ab);
    MAP_inert(response+1) = integral2(fun_MAP_InertPeak, t_min, t_max, s_min, s_max)/MAP_norm;
end

% Expected inertial (to be minimized):
MAP_px = normcdf(next_x, t_estim, s_estim);
MAP_exp_inert = (1 - MAP_px) * MAP_inert(1) + MAP_px * MAP_inert(2);

end
