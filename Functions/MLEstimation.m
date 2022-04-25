function result = MLEstimation(params, version)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Maximum likelihood estimation methods (Modular, probability Space)
%   Different versions for selecting next stimulation amplitude
%   Version 1: Awiszus, 2003 
%   Version 2: Full MLE explicit
%   Version 3: MLE
%	Version 4: MAP
%	Version 5: MAP min variance
%	Version 6: MAP min inertia
%   Version 7: Awiszus, 2003, v2
%   Version 8: Full MLE explicit v2

if nargin == 1
    version = 1;
end

step_number = params.step_number;
amplitude_list = Inf(1, step_number);
response_list = Inf(1, step_number);
response_bin = false(1, step_number);

MLE_t = NaN(1, step_number);
MLE_s = NaN(1, step_number);
MLE_t_mean = zeros(1, step_number);
MLE_t_var = zeros(1, step_number);
MAP_t = NaN(1, step_number);
MAP_s = NaN(1, step_number);
MAP_t_mean = zeros(1, step_number);
MAP_t_var = zeros(1, step_number);

amplitude_list(1) = params.start_amplitude;

MLE_opts = params.opts;
MAP_opts = params.opts;

t_min = 0.01;
t_max = 1.2;
s_min = 0.001;
s_max = 0.5;
       
switch version
    case 1
        amplitude_step = 0.001;
        possible_amplitude_list = 0.20 : amplitude_step : 1.05;
        
    case 3
        amplitude_step = 0.001;
        possible_amplitude_list = 0.20 : amplitude_step : 1.05;
        MLE_grid_size = 401;
        possible_slope_list = logspace(log10(s_min), log10(s_max), MLE_grid_size);
        [possible_amplitude_grid, possible_slope_grid] = ndgrid(possible_amplitude_list, possible_slope_list);
end

             
for step_cnt = 1 : step_number
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Stimulate with next (chosen) amplitude
    response_list(step_cnt) = virtstimulate(amplitude_list(step_cnt), params.subj_parameters);
    response_bin(step_cnt) = logical(response_list(step_cnt) > params.y_thresh);
    
    stim_vec = [0, amplitude_list(1:step_cnt)]; % added 0-stimulation (with 0 response)
    response_vec = [false, response_bin(1:step_cnt)];
    
    switch version
        case {1,7}
            LB = t_min;   % lower bounds for threshold in MLE
            UB = t_max;      % upper bounds for threshold in MLE

            %%%%%%%%%%%%%%%%%%%%%%
            %   MLE
            fun_MLE_min = @(theta) -1 * loglikelyhood(stim_vec, response_vec, theta, theta*0.07);
            MLE_init = 0.5;
            %[MLE, ~, MLE_exitflag] = fminsearch(fun_MLE_min, MLE_init, MLE_opts);	% find max via min of negative
            [MLE, ~, MLE_exitflag] = fmincon(fun_MLE_min, MLE_init, [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
            if ~MLE_exitflag            % Optimization did not converge, do it again:
                count_run = 1;
                %fprintf('MLE not successful. Rerun No. %d...\n', count_run);
                MLE_init = rand * 1;
                [MLE, ~, MLE_exitflag] = fmincon(fun_MLE_min, MLE_init, [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
                while ~MLE_exitflag        % do it again:
                    count_run = count_run + 1;
                    %fprintf('MLE not successful again. Rerun No. %d...\n', count_run);
                    MLE_init = rand * 1;
                    [MLE, ~, MLE_exitflag] = fmincon(fun_MLE_min, MLE_init, [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
                end
            end
            
            if MLE_exitflag
                MLE_t(step_cnt) = MLE;
                MLE_s(step_cnt) = MLE * 0.07;
            end
            
            %   Error/Deviation;  Termination Information
            %   Normalization of probability density:
            fun_MLE_E0 = @(t) pMLE(stim_vec, response_vec, t, t*0.07);
            p_normalization = integral(fun_MLE_E0, t_min, t_max);
            %   1st order momentum:
            fun_MLE_E1 = @(t) t .* pMLE(stim_vec, response_vec, t, t*0.07);
            MLE_t_mean(step_cnt) = integral(fun_MLE_E1, t_min, t_max) / p_normalization;
            %   2nd order momentum and variance:
            % MLE_E2 = @(t,s) (t.^2).*pMLE(stim_vec, response_vec, t, s);
            fun_MLE_var = @(t) ((t-MLE_t_mean(step_cnt)).^2).*pMLE(stim_vec, response_vec, t, t*0.07);
            MLE_t_var(step_cnt) = integral(fun_MLE_var, t_min, t_max) / p_normalization;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %   MAP
            fun_MAP_min = @(theta) (-1) * logMAP(stim_vec, response_vec, theta, theta*0.07, @p_ab);
            MAP_init = 0.5;
            %[MAP, ~, MAP_exitflag] = fminsearch(fun_MAP_min, MAP_init, MAP_opts);	% find max via min of negative
            [MAP, ~, MAP_exitflag] = fmincon(fun_MAP_min, MAP_init, [], [], [], [], LB, UB, [], MAP_opts);  % find max via min of negative
            if ~MAP_exitflag
                count_run = 1;
                %fprintf('Step %d: MAP not successful. Rerun No. %d...\n', step_cnt, count_run);
                MAP_init = rand * 1;
                [MAP, ~, MAP_exitflag] = fmincon(fun_MAP_min, MAP_init, [], [], [], [], LB, UB, [], MAP_opts);  % find max via min of negative
                while ~MAP_exitflag
                    count_run = count_run + 1;
                    %fprintf('Step %d: MAP not successful again. Rerun No. %d...\n', step_cnt, count_run);
                    MAP_init = rand * 1;
                    [MAP, ~, MAP_exitflag] = fmincon(fun_MAP_min, MAP_init, [], [], [], [], LB, UB, [], MAP_opts);  % find max via min of negative
                end
            end
            if MAP_exitflag
                MAP_t(step_cnt) = MAP;
                MAP_s(step_cnt) = MAP * 0.07;
            end
            
            %   Error/Deviation;  Termination Information
            %   Normalization of probability density:
            fun_MAP_E0 = @(t) pMAP(stim_vec, response_vec, t, t*0.07, @p_ab);
            p_normalization = integral(fun_MAP_E0, t_min, t_max);
            %   1st order momentum:
            fun_MAP_E1 = @(t) t.*pMAP(stim_vec, response_vec, t, t*0.07, @p_ab);
            MAP_t_mean(step_cnt) = integral(fun_MAP_E1, t_min, t_max) / p_normalization;
            %   2nd order momentum and variance:
            % MAP_E2 = @(t,s) (t.^2).*pMAP(stim_vec, response_vec, t, s, @myp_ab_fct);
            fun_MAP_var = @(t) ((t-MAP_t_mean(step_cnt)).^2).*pMAP(stim_vec, response_vec, t, t*0.07, @p_ab);
            MAP_t_var(step_cnt) = integral(fun_MAP_var, t_min, t_max) / p_normalization;
            
        otherwise
            %%%%%%%%%%%%%%%%%%%%%%
            %   MLE
            LB = [t_min, s_min];   % lower bounds for threshold and slope in MLE
            UB = [t_max, s_max];      % upper bounds for threshold and slope in MLE

            fun_MLE_min = @(theta) -1 * loglikelyhood(stim_vec, response_vec, theta(1), theta(2));
            MLE_init = [0.5, 0.2];
            %[MLE, ~, MLE_exitflag] = fminsearch(fun_MLE_min, MLE_init, MLE_opts);	% find max via min of negative
            [MLE, ~, MLE_exitflag] = fmincon(fun_MLE_min, MLE_init, [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
            if ~MLE_exitflag            % Optimization did not converge, do it again:
                count_run = 1;
                %fprintf('MLE not successful. Rerun No. %d...\n', count_run);
                MLE_init = [rand * 1, rand * 0.15];
                [MLE, ~, MLE_exitflag] = fmincon(fun_MLE_min, MLE_init, [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
                while ~MLE_exitflag        % do it again:
                    count_run = count_run + 1;
                    %fprintf('MLE not successful again. Rerun No. %d...\n', count_run);
                    MLE_init = [rand * 1, rand * 0.1];
                    [MLE, ~, MLE_exitflag] = fmincon(fun_MLE_min, MLE_init, [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
                end
            end
            
            if MLE_exitflag
                MLE_t(step_cnt) = MLE(1);
                MLE_s(step_cnt) = MLE(2);
            end
            
            %   Error/Deviation;  Termination Information
            %   Normalization of probability density:
            fun_MLE_E0 = @(t,s) pMLE(stim_vec, response_vec, t, s);
            p_normalization = integral2(fun_MLE_E0, t_min, t_max, s_min, s_max);
            %   1st order momentum:
            fun_MLE_E1 = @(t,s) t .* pMLE(stim_vec, response_vec, t, s);
            MLE_t_mean(step_cnt) = integral2(fun_MLE_E1, t_min, t_max, s_min, s_max) / p_normalization;
            %   2nd order momentum and variance:
            % MLE_E2 = @(t,s) (t.^2).*pMLE(stim_vec, response_vec, t, s);
            fun_MLE_var = @(t,s) ((t-MLE_t_mean(step_cnt)).^2).*pMLE(stim_vec, response_vec, t, s);
            MLE_t_var(step_cnt) = integral2(fun_MLE_var, t_min, t_max, s_min, s_max) / p_normalization;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %   MAP
            fun_MAP_min = @(theta) (-1) * logMAP(stim_vec, response_vec, theta(1), theta(2), @p_ab);
            MAP_init = [0.5, 0.15];
            %[MAP, ~, MAP_exitflag] = fminsearch(fun_MAP_min, MAP_init, MAP_opts);	% find max via min of negative
            [MAP, ~, MAP_exitflag] = fmincon(fun_MAP_min, MAP_init, [], [], [], [], LB, UB, [], MAP_opts);  % find max via min of negative
            if ~MAP_exitflag
                count_run = 1;
                %fprintf('Step %d: MAP not successful. Rerun No. %d...\n', step_cnt, count_run);
                MAP_init = [rand * 1, rand * 0.10];
                [MAP, ~, MAP_exitflag] = fmincon(fun_MAP_min, MAP_init, [], [], [], [], LB, UB, [], MAP_opts);  % find max via min of negative
                while ~MAP_exitflag
                    count_run = count_run + 1;
                    %fprintf('Step %d: MAP not successful again. Rerun No. %d...\n', step_cnt, count_run);
                    MAP_init = [rand * 1, rand * 0.05];
                    [MAP, ~, MAP_exitflag] = fmincon(fun_MAP_min, MAP_init, [], [], [], [], LB, UB, [], MAP_opts);  % find max via min of negative
                end
            end
            if MAP_exitflag
                MAP_t(step_cnt) = MAP(1);
                MAP_s(step_cnt) = MAP(2);
            end
            
            %   Error/Deviation;  Termination Information
            %   Normalization of probability density:
            fun_MAP_E0 = @(t,s) pMAP(stim_vec, response_vec, t, s, @p_ab);
            p_normalization = integral2(fun_MAP_E0, t_min, t_max, s_min, s_max);
            %   1st order momentum:
            fun_MAP_E1 = @(t,s) t.*pMAP(stim_vec, response_vec, t, s, @p_ab);
            MAP_t_mean(step_cnt) = integral2(fun_MAP_E1, t_min, t_max, s_min, s_max) / p_normalization;
            %   2nd order momentum and variance:
            % MAP_E2 = @(t,s) (t.^2).*pMAP(stim_vec, response_vec, t, s, @myp_ab_fct);
            fun_MAP_var = @(t,s) ((t-MAP_t_mean(step_cnt)).^2).*pMAP(stim_vec, response_vec, t, s, @p_ab);
            MAP_t_var(step_cnt) = integral2(fun_MAP_var, t_min, t_max, s_min, s_max) / p_normalization;
    end
    
    
    if (step_cnt < step_number)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Next Amplitude
        switch version
            case 1
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %   Awiszus with 0.07 parameter
                fun_awi_FindNextAmplMeasure = @(t) (-1)*loglikelyhood(stim_vec, response_vec, t, 0.07*t); 
                
                possible_amplitude_likelihood = fun_awi_FindNextAmplMeasure(possible_amplitude_list);
                ampl_inds = possible_amplitude_likelihood <= min(possible_amplitude_likelihood);
                next_amplitude = mean(possible_amplitude_list( ampl_inds ));		
                % minimum can be a whole list of values, if several thresholds have same likelihood => take mean
            case 2
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %   Full MLE explicit, with mean if more than 1 solution:
                %   Place the grid around the minimum found and only evaluate there, if more than one:
                fun_MLE_FindNextAmplMeasure = @(t,s) (-1)*loglikelyhood(stim_vec, response_vec, t, s);
                
                possible_amplitude_likelihoodfull = fun_MLE_FindNextAmplMeasure(possible_amplitude_grid, possible_slope_grid);
                [ampl_inds, ~] = find( possible_amplitude_likelihoodfull <= min(possible_amplitude_likelihoodfull, [],  'all'));
                next_amplitude = mean(possible_amplitude_list( ampl_inds ));
                % minimum can be a whole list of values, if several thresholds have same likelihood => take mean
            case 3
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %   MLE
                next_amplitude = MLE_t(step_cnt);
            case 4
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %   MAP
                next_amplitude = MAP_t(step_cnt);
            case 5
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %   MAP Minimize Variance of parameter t
                fun_MAP_var_min = @(next_x) MAPexpectedvar(stim_vec, response_vec, next_x, MAP_t(step_cnt), MAP_s(step_cnt), @p_ab);
                % next_amplitude = fminsearch(fun_MAP_var_min, MAP_t(step_cnt), MAP_opts);
                next_amplitude = fmincon(fun_MAP_var_min, MAP_t(step_cnt), [], [], [], [], t_min, t_max, [], MAP_opts);  % find max via min of negative
            case 6
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %   MAP Minimize Inertia of parameter t around its peak
                fun_MAP_inert_min = @(next_x) MAPexpectedinertpeak(stim_vec, response_vec, next_x, MAP_t(step_cnt), MAP_s(step_cnt), @p_ab);
                % next_amplitude = fminsearch(fun_MAP_inert_min, MAP_t(step_cnt), MAP_opts);
                next_amplitude = fmincon(fun_MAP_inert_min, MAP_t(step_cnt), [], [], [], [], t_min, t_max, [], MAP_opts);  % find max via min of negative
            case 7
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %   Awiszus with 0.07 parameter v2, no averaging
                fun_awi_FindNextAmplMeasure = @(m) (-1)*loglikelyhood(stim_vec, response_vec, m, 0.07*m); 
                % next_amplitude= fminsearch(fun_awi_FindNextAmplMeasure, MLE_t(step_cnt), MLE_opts);        % find max via min of negative
                next_amplitude = fmincon(fun_awi_FindNextAmplMeasure, MLE_t(step_cnt),  [], [], [], [], t_min, t_max, [], MAP_opts);  % find max via min of negative
            case 8
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %   Full MLE explicit v2, no averaging
                fun_MLE_FindNextAmplMeasure = @(theta) (-1)*loglikelyhood(stim_vec, response_vec, theta(1), theta(2));
                % tmp = fminsearch(fun_MLE_FindNextAmplMeasure, [MLE_t(step_cnt), MLE_s(step_cnt)], MLE_opts);        % find max via min of negative
                tmp = fmincon(fun_MLE_FindNextAmplMeasure, [MLE_t(step_cnt), MLE_s(step_cnt)], [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
                next_amplitude = tmp(1);
        end

        amplitude_list(step_cnt+1) = max(next_amplitude, 0);
    end

end

MLE_abs_err = MLE_t - params.thresh_x;              % absolute error (0-1)
MLE_rel_err = MLE_abs_err ./ params.thresh_x *100;  % relative error (0%-100%)

MAP_abs_err = MAP_t - params.thresh_x;             	% absolute error (0-1)
MAP_rel_err = MAP_abs_err ./ params.thresh_x *100;  % relative error (0%-100%)


result = struct('amplitude_list', amplitude_list, ...
                'response_list', response_list, ...
                'MLE_t', MLE_t, ...
                'MLE_s', MLE_s, ...
                'MLE_t_mean', MLE_t_mean, ...
                'MLE_t_var', MLE_t_var ,...
                'MLE_abs_err', MLE_abs_err,...
                'MLE_rel_err', MLE_rel_err,...
                'MAP_t', MAP_t, ...
                'MAP_s', MAP_s, ...
                'MAP_t_mean', MAP_t_mean, ...
                'MAP_t_var', MAP_t_var,...
                'MAP_abs_err', MAP_abs_err,...
                'MAP_rel_err', MAP_rel_err...
               );
end
