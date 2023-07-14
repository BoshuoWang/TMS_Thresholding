function result = MLEstimation(params, version)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Maximum likelihood estimation methods (Modular, probability space)
%   Different versions, rewritten Nov. 2022
%   Version 1: MLE-1 estimation & stepping, explicit evaluation (Awiszus, 2003) 
%   Version 2: MLE-1 estimation & stepping, direct maximization
%   Version 3: MAP-1 estimation & stepping, direct maximization
%
%   Version 4: MLE-2 estimation & stepping, explicit evaluation
%	Version 5: MLE-2 estimation & stepping, direct maximization
%   Version 6: MAP-2 estimation & stepping, direct maximization
%
%	Version 7: MAP-2 estimation, direct maximization, min variance stepping
%	Version 8: MAP-2 estimation, direct maximization, min inertia stepping  

if nargin == 1
    version = 1;
end

MLE_opts = params.opts;
MAP_opts = params.opts;

t_step = 0.002;
t_min = t_step;     t_max = 1.3;        % Threshold parameter can be larger than 100% MSO
s_min = 0.005;      s_max = 0.5;        % Two orders-of-magnitude

step_number = params.step_number;
amplitude_list = NaN(1, step_number + 2);   
response_list  = NaN(1, step_number + 2);
response_bin = false(1, step_number + 2);
thresh_est     = NaN(1, step_number + 2);
run_time       = NaN(1, step_number + 2);
% Add default samples (pseudo-responses) at end of range: 
% no response at 0% MSO, and response at t_max
amplitude_list(1) = 0;              response_bin(1) = false;
amplitude_list(2) = t_max + t_step; response_bin(2) = true;
amplitude_list(3) = min(params.start_amplitude, t_max);



if version <= 3                 % Single parameter MLE/MAP
    LB = t_min;                 % Lower bounds for threshold 
    UB = t_max;                 % Upper bounds for threshold
else                            % Two-parameter MLE/MAP
    LB = [t_min, s_min];        % Lower bounds for threshold and spread
    UB = [t_max, s_max];        % Upper bounds for threshold and spread
end

if version == 1                 % Amplitude list for explicit evaluation
    possible_amplitude_list = t_min : t_step : t_max;
end
if version == 4                 % Grid of amplitudes and spreads for explicit evaluation
    possible_amplitude_list = t_min : t_step : t_max;
    MLE_grid_size = round(log10(s_max) - log10(s_min)) * 50;     % 50 points per decade, 100 points total
    tmp_vec = logspace(log10(s_min), log10(s_max), MLE_grid_size*2+1);	
    possible_spread_list = tmp_vec(2:2:end-1);         % Sample at center of intervals
    [possible_amplitude_grid, possible_spread_grid] = ndgrid(possible_amplitude_list, possible_spread_list);
end

for step_cnt = (1 : step_number) + 2    % offset by 2 to include default samples
    T_start = tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Stimulate with next (chosen) amplitude
    response_list(step_cnt) = virtstimulate(amplitude_list(step_cnt), params.subj_parameters);
    response_bin(step_cnt)  = logical(response_list(step_cnt) > params.y_thresh);
    
    vec_ind = 1 : step_cnt;             % include default samples
    
    switch version                      % function to maximize
        case {1,2}                      % Versions 1,2: MLE-1
            fun_min = @(t)      (-1)*loglikelyhood(amplitude_list(vec_ind), response_bin(vec_ind), t,        0.07*t);
        case {3}                        % Version 3: MAP-1
            fun_min = @(t)      (-1)*logMAP(       amplitude_list(vec_ind), response_bin(vec_ind), t,        0.07*t,   @p_ab);
        case {4}                        % Version 4: MLE-2, explicit
            fun_min = @(t,s)    (-1)*loglikelyhood(amplitude_list(vec_ind), response_bin(vec_ind), t,        s);
        case {5}                        % Version 5: MLE-2
            fun_min = @(theta)  (-1)*loglikelyhood(amplitude_list(vec_ind), response_bin(vec_ind), theta(1), theta(2));
        case {6,7,8}                    % Versions 6, 7, 8: MAP-2
            fun_min = @(theta)  (-1)*logMAP(       amplitude_list(vec_ind), response_bin(vec_ind), theta(1), theta(2), @p_ab);
    end    
    
    switch version                      % Maximization step for threshold estimation
        case {1}                        % Version 1: explicit likelihood evaluation, 1 parameter
            explicit_likelihood = fun_min(possible_amplitude_list);
            ampl_inds = (explicit_likelihood <= min(explicit_likelihood));          % maxima can be a whole list of values
            thresh_est(step_cnt) = mean(possible_amplitude_list( ampl_inds ));      % if several thresholds have same likelihood => take mean
        case {2,3}                      % Versions 2,3: direct maximization, 1 parameter
            count_run = 1;
            theta_init = amplitude_list(step_cnt);
            [theta_min, ~, exitflag] = fmincon(fun_min, theta_init,  [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
            % [X,FVAL,EXITFLAG] = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
            while exitflag <= 0            % Optimization did not converge, do it again:
                fprintf('MLE-1/MAP-1 not successful on run No. %d. Rerunning...\n', count_run);
                count_run = count_run + 1;
                theta_init = rand * t_max;
                [theta_min, ~, exitflag] = fmincon(fun_min, theta_init,  [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
            end
            thresh_est(step_cnt) = theta_min;
        case {4}                        % Version 4: explicit likelihood evaluation, 2 parameters
            explicit_likelihood = fun_min(possible_amplitude_grid, possible_spread_grid);
            [ampl_inds, ~] = find( explicit_likelihood <= min(explicit_likelihood, [],  'all'));    % maxima can be a whole list of values
            thresh_est(step_cnt) = mean(possible_amplitude_list( ampl_inds ));                      % if several thresholds have same likelihood => take mean
        case {5,6,7,8}                  % Version 5,6,7,8: direct maximization, 2 parameters
            count_run = 1;
            theta_init = [amplitude_list(step_cnt), sqrt(s_min*s_max)];
            [theta_min, ~, exitflag] = fmincon(fun_min, theta_init, [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
            % [X,FVAL,EXITFLAG] = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
            while exitflag <= 0            % Optimization did not converge, do it again:
                fprintf('MLE-2/MAP-2 not successful on run No. %d. Rerunning...\n', count_run);
                count_run = count_run + 1;
                theta_init = [rand * t_max, rand * s_max];
                [theta_min, ~, exitflag] = fmincon(fun_min, theta_init,  [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
            end
            thresh_est(step_cnt) = theta_min(1);
    end
   
    if version <= 6                 % Stepping rule: versions 1-6, use threshold
        next_amplitude = thresh_est(step_cnt);
    else                            % Stepping rule: versions 7,8, minimize expected second-order function
        if version == 7
            fun_MAP_next_min = @(next_x) MAPexpectedvar(      amplitude_list(vec_ind), response_bin(vec_ind), next_x, theta_min(1), theta_min(2), @p_ab);
        else
            fun_MAP_next_min = @(next_x) MAPexpectedinertpeak(amplitude_list(vec_ind), response_bin(vec_ind), next_x, theta_min(1), theta_min(2), @p_ab);
        end
        [next_amplitude, ~, exitflag] = fmincon(fun_MAP_next_min, theta_min(1), [], [], [], [], t_min, t_max, [], MAP_opts);  % find max via min of negative
        while exitflag <= 0
            fprintf('Fmincon not successful on run No. %d. Rerunning...\n', count_run);
            [next_amplitude, ~, exitflag] = fmincon(fun_MAP_next_min, rand * t_max, [], [], [], [], t_min, t_max, [], MAP_opts);  % find max via min of negative
        end
    end
    if (step_cnt < step_number + 2)
        amplitude_list(step_cnt+1) = next_amplitude;    % Stimulation amplitude automatically clipped between t_min, t_max (0% and 130% MSO)
    end
    run_time(step_cnt) = toc(T_start);
end

abs_err = thresh_est - params.thresh_x;         % absolute error
rel_err = abs_err ./ params.thresh_x *100;      % relative error 

result = struct('amplitude_list', amplitude_list(3:end), ...
                'response_list', response_list(3:end), ...
                'thresh_est', thresh_est(3:end), ...
                'abs_err', abs_err(3:end), ...
                'rel_err', rel_err(3:end), ...
                'run_time', run_time(3:end) ...
               );                               % remove default samples and cut to right size
end
