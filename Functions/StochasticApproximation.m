%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   StochasticApproximation, Robbins-Monro, fully log, all versions
%   Version 1, 2:       non-adaptive 1/i, 1st & 2nd order (original versions 9, 10)
%   Version 3, 4:       adaptive 1/i, 1st & 2nd order (original versions 1/11/81, 2)
%   Version 5, 6:       adaptive 2^(-i), 1st & 2nd order (original versions 3, 4)
%   Version 7, 8:       adaptive 2^(-i/2), 1st & 2nd order (original versions 5/15, 6)
%   Version 9, 10:      adaptive 2^(-i/2), 1st & 2nd order, de- and increasing steps (original versions 7, 8)
%   Version 11, 12:     adaptive 1/i, 1st & 2nd order, step increases if plateau reached (original version 21)
%   Version 13:         Stochastic Newton, starts as RM (adaptive 1/i, 1st order) until robust estimation with linear regression (orginal version snewta)
%   Version 7-3:        adaptive 2^(-i/2) switches to 1/i, 1st order i.e., rapid approach switching to a.s. convergence. Transition: max step stize < +/- 0.015 = 1.5% MSO (original version 151)

%  StochasticApproximation(params, version, is_analog, num_start_conditions, num_second_order_weights, run_lin, run_MLE)
%  Default argin after params: 1, true, 1, 1, false, false

%   Estimate derivative with robust linear regression
%   Run MLE in parallel

function result = StochasticApproximation(...
    params, version, is_analog, ...
    num_start_conditions, num_second_order_weights, ...
    run_lin, run_MLE)

if nargin < 2
    version =  1;
end
if nargin < 3
    is_analog = true;
end
if nargin < 4
    num_start_conditions = 1;
end
% if is_analog
    default_start_ctrl_seqs = 0.20/(log10(1e-2)-log10(1e-5));   
% else
%     default_start_ctrl_seqs = 0.20/2;   
% end
if num_start_conditions == 1
    start_ctrl_seqs = default_start_ctrl_seqs;   
    % starting value is approx. cotangent (ai*dy = dx => ai = dx/dy)
    % Overestimate slightly??
else
    num_start_conditions = fix(num_start_conditions/2)*2 + 1; % Make odd 
    start_ctrl_seqs = default_start_ctrl_seqs * logspace(-1, 1, num_start_conditions);
end
if nargin < 5
    num_second_order_weights = 1;
end
if num_second_order_weights == 1
    second_order_weights = -0.10;
else
    num_second_order_weights = fix(num_second_order_weights/2)*2 + 1; % Make odd 
    second_order_weights = linspace(-1, 1, num_second_order_weights);
end
num_conditions = num_start_conditions * num_second_order_weights;

if nargin < 6 || (num_conditions > 1)
    run_lin = false;
end
if nargin < 7 || (num_conditions > 1)
    run_MLE = false;
end


y_thresh = params.y_thresh;
step_number = params.step_number;

[start_ctrl_seqs, second_order_weights] = ndgrid(start_ctrl_seqs, second_order_weights);

control_sequence = NaN(step_number, num_conditions);
amplitude_list = NaN(step_number+1, num_conditions);
response_bin = false(step_number, num_conditions);
if is_analog
    response_list = NaN(step_number, num_conditions);
end

% delta_ctrl_seq = zeros(1, num_conditions);
delta_Y = zeros(step_number, num_conditions);

stochastic_mode = true(step_number, 1);

amplitude_list(1, :) = params.start_amplitude;

% (compensating that saturation reduces the error term and thus the drive to go back to the threshold)
cutoff_low = 10e-6;            % 10 µV. below that, the step size is increased again 
cutoff_high = 10e-3;           % 10 mV. above that, the step size is increased again

if run_lin
    lin_estim = NaN(1, step_number);
    min_lin = 5;
end
if run_MLE
    MLE_t = NaN(1, step_number);
    MLE_s = NaN(1, step_number);
    MLE_opts = params.opts;
end

switch version
    case {3,4, 11,12, 13,14}
        number_of_sign_changes = zeros(1, num_conditions);
    case {5,6}
        b_i = 2 ;
    case {7,8, 9,10}
        b_i = sqrt(2);
    case {73,84}
        b_i = sqrt(2);
        as_convergence = false(1, num_conditions);   % false: geometric; true: harmonic
        number_of_sign_changes = zeros(1, num_conditions);
end


for step_cnt = 1 : step_number
    response =  virtstimulate(amplitude_list(step_cnt, :), params.subj_parameters) ;
    % response > 0; real function is already taken in virtstimulate
    response_bin(step_cnt, :) = logical(response >= y_thresh);
    if is_analog
        response_list(step_cnt, :) = response;
        delta_Y(step_cnt, :) = log10(response) - log10(y_thresh); % in log-space
    else
        delta_Y(step_cnt, :) = 2 * response_bin(step_cnt, :) - 1; % convert from logical [0,1] to [-1,1]
    end
    if run_lin
        %%%%%%%%%%%%%%%%%%%%%%
        % Estimate Derivative:
        if (step_cnt > min_lin)
            [pcoeff, ~] = robustfit(amplitude_list(1:step_cnt), log10(response_list(1:step_cnt)), 'bisquare');
            % Estimate threshold from this fit:
            lin_estim(step_cnt) = (log10(y_thresh) - pcoeff(1)) / pcoeff(2);
            % a + b*x = y -> x = (y-a)/b
        end
    end
    if run_MLE
        %%%%%%%%%%%%%%%%%%%%%%
        % MLEstimator with data
        stim_vec = cat(2, zeros(1, num_conditions)', amplitude_list(1:step_cnt, :)'); % added 0-stimulation (with 0 response)
        response_vec = cat(2, false(1, num_conditions)', response_bin(1:step_cnt, :)');

        fun_MLE_min = @(theta) -1 * loglikelyhood(stim_vec, response_vec, theta(1), theta(2));
        MLE_init = [0.5, 0.2];
        [MLE, ~, MLE_exitflag] = fminsearch(fun_MLE_min, MLE_init, MLE_opts);	% find max via min of negative
        if ~MLE_exitflag            % Optimization did not converge, do it again:
            MLE_init = [rand *1, rand *0.15];
            [MLE, ~, MLE_exitflag] = fminsearch(fun_MLE_min, MLE_init, MLE_opts);
            while ~MLE_exitflag        % do it again:
                MLE_init = [rand *1, rand *0.1];
                [MLE, ~, MLE_exitflag] = fminsearch(fun_MLE_min, MLE_init, MLE_opts);
            end
        end
        if MLE_exitflag
            MLE_t(step_cnt) = MLE(1);
            MLE_s(step_cnt) = MLE(2);
        end
    end
    
    % Calculate control sequence and next amplitude
    switch version
        case {1,2}      % ~1/i, original version IX, X
            control_sequence(step_cnt, :) = start_ctrl_seqs(:)' / step_cnt;
        case {3,4}      % overall: ~1/i, only in case of sign change, original version I, II
            if step_cnt == 1
                control_sequence(step_cnt, :) = start_ctrl_seqs(:)';
            else        % step_cnt > 1
                sign_change = xor(response_bin(step_cnt, :), response_bin(step_cnt-1, :));  % logical, 0: no sign change, 1: sign change
                control_sequence(step_cnt, ~sign_change) = control_sequence(step_cnt-1, ~sign_change);
                
                number_of_sign_changes = number_of_sign_changes + sign_change;              % none-zero for those with sign changes
                control_sequence(step_cnt,  sign_change) = control_sequence(step_cnt-1,  sign_change) .* number_of_sign_changes(sign_change)./ (number_of_sign_changes(sign_change) + 1);
            end
        case {5,6,7,8}      % overall: ~2^(-i) or ~2^(-i/2), only in case of sign change, original version III, IV
            if step_cnt == 1
                control_sequence(step_cnt, :) = start_ctrl_seqs(:)';
            else % step_cnt > 1
                sign_change = xor(response_bin(step_cnt, :), response_bin(step_cnt-1, :));      % logical, 0: no sign change, 1: sign change
                control_sequence(step_cnt, :) = control_sequence(step_cnt-1, :) ./ (b_i.^sign_change);  %  b_i^-1 for sign change
            end
        case {9,10}     % overall: ~2^(-i/2), only in case of sign change, original version VII, VIII
%             delta_ctrl_seq = (b_i    - 1) .*  xor(response_bin(step_cnt, :), response_bin(step_cnt-1, :)) + ...
%                              (b_i^-1 - 1) .* ~xor(response_bin(step_cnt, :), response_bin(step_cnt-1, :)) ;
            if step_cnt == 1
                control_sequence(step_cnt, :) = start_ctrl_seqs(:)';
            else % step_cnt > 1
                sign_change = xor(response_bin(step_cnt, :), response_bin(step_cnt-1, :)) * 2 -1;   % -1: no sign change, 1: sign change
                control_sequence(step_cnt, :) = control_sequence(step_cnt-1, :) ./ (b_i.^sign_change);
            end
        case {11,12}    % overall: ~1/i, only in case of sign change & not in plateau regions, original version XXI
            if step_cnt == 1
                control_sequence(step_cnt, :) = start_ctrl_seqs(:)';
            else        % step_cnt > 1
                sign_change = xor(response_bin(step_cnt, :), response_bin(step_cnt-1, :)) & ...
                                ( response_list(step_cnt, :) > cutoff_low  )  & ...
                                ( response_list(step_cnt, :) < cutoff_high );       % logical, 1: sign change, counted only outside plateau, 0: no sign change or sign change but in plateau
                control_sequence(step_cnt, ~sign_change) = control_sequence(step_cnt-1, ~sign_change);
                
                number_of_sign_changes = number_of_sign_changes + sign_change;                     
                control_sequence(step_cnt,  sign_change) = control_sequence(step_cnt-1,  sign_change) .* number_of_sign_changes(sign_change)./ (number_of_sign_changes(sign_change) + 1);
            end
        case {13,14} % Stochastic Newton, starts as RM (adaptive 1/n, 1st order) until robust estimation with linear regression, original version swneta
            % Linear approximation considered stable as soon as
            % derivative positive and steeper than initial derivative
            % as well as estimated threshold between 10% and 100%
            % amplitude => Newton-Raphson-like descent using the
            % linearized approximation
            if (step_cnt > min_lin) && (pcoeff(2) > 1/start_ctrl_seqs) && ...
               (lin_estim(step_cnt)>0.1) && (lin_estim(step_cnt)<1.0)
                stochastic_mode(step_cnt) = false;
            else    % Normal adaptive RM until robust estimate
                if step_cnt == 1
                    control_sequence(step_cnt, :) = start_ctrl_seqs(:)';
                else        % step_cnt > 1
                    sign_change = xor(response_bin(step_cnt, :), response_bin(step_cnt-1, :));  % logical, 0: no sign change, 1: sign change
                    control_sequence(step_cnt, ~sign_change) = control_sequence(step_cnt-1, ~sign_change);
                    
                    number_of_sign_changes = number_of_sign_changes + sign_change;              % none-zero for those with sign changes
                    control_sequence(step_cnt,  sign_change) = control_sequence(step_cnt-1,  sign_change) .* number_of_sign_changes(sign_change)./ (number_of_sign_changes(sign_change) + 1);
                end
            end
            
        case {73,84}	% Controlling sequence in the beginning: sqrt(1/2)^n; only in case of sign change)
            % switches to 1/n as soon as max expected step size < +/- 1.5% MSO
            if step_cnt == 1
                control_sequence(step_cnt, :) = start_ctrl_seqs(:)';
            else        % step_cnt > 1
                sign_change = xor(response_bin(step_cnt, :), response_bin(step_cnt-1, :));  % logical, 0: no sign change, 1: sign change
                control_sequence(step_cnt, ~sign_change) = control_sequence(step_cnt-1, ~sign_change);
                
                max_step_size = control_sequence(step_cnt-1, :) .* ...
								(    log10(max(response_list, [], 1, 'omitnan') ) - ...
								 min(log10(min(response_list, [], 1, 'omitnan') ), log10(y_thresh)) );
                as_convergence = as_convergence | (max_step_size <  1.5e-2);            % switch mode and stick to a.s.; false: geometric; true: harmonic
                
                ind_sign_change_NOT_as = sign_change & ~as_convergence;
                control_sequence(step_cnt,  ind_sign_change_NOT_as) = control_sequence(step_cnt-1,  ind_sign_change_NOT_as) / (b_i);
                
                ind_sign_change_AND_as = sign_change &  as_convergence;
                number_of_sign_changes = number_of_sign_changes + sign_change;              % none-zero for those with sign changes
                control_sequence(step_cnt,  ind_sign_change_AND_as) = control_sequence(step_cnt-1,  ind_sign_change_AND_as) .* number_of_sign_changes(ind_sign_change_AND_as)./ (number_of_sign_changes(ind_sign_change_AND_as) + 1);
            end
    end
    
    % Finding next amplitude
    if stochastic_mode(step_cnt)
        % x_(i+1) = x_i - a_i * (y_i - y_th)
        amplitude_list(step_cnt+1, :) = amplitude_list(step_cnt, :) - ...
                                        control_sequence(step_cnt, :) .* delta_Y(step_cnt, :);
        if ~mod(version, 2) && (step_cnt >= 2)   % Even versions, add second order term
            % x_(i+1) = x_i - a_i * ( (y_i - y_th) + w * (y_(i-1) - y_th) )
            amplitude_list(step_cnt+1, :) = amplitude_list(step_cnt+1, :) - ...
                                            control_sequence(step_cnt, :) .* (second_order_weights(:)' .* delta_Y(step_cnt-1, :) ) ; 
        end
        amplitude_list(step_cnt+1, :) = max(0, amplitude_list(step_cnt+1, :)); % no negative stimuli
    else
        amplitude_list(step_cnt+1) = lin_estim(step_cnt);
    end
end


result.control_sequence = squeeze(reshape(control_sequence, [step_number, num_start_conditions, num_second_order_weights]));
result.amplitude_list = squeeze(reshape(amplitude_list, [step_number+1, num_start_conditions, num_second_order_weights]));
result.step_size = squeeze(reshape(cat(1, zeros(1, num_conditions), diff(amplitude_list)), [step_number+1, num_start_conditions, num_second_order_weights]));

result.abs_err = result.amplitude_list - params.thresh_x;               % absolute error (0-1)
result.rel_err = result.abs_err ./ params.thresh_x *100;                % relative error (0%-100%)

if is_analog
    result.response_list = squeeze(reshape(response_list, [step_number, num_start_conditions, num_second_order_weights]));
else
    result.response_bin = squeeze(reshape(response_bin, [step_number, num_start_conditions, num_second_order_weights]));
end
if run_lin
    result.lin_estim = squeeze(reshape(lin_estim, [step_number, num_start_conditions, num_second_order_weights]));
    result.lin_abs_err = result.lin_estim - params.thresh_x;            % absolute error (0-1)
    result.lin_rel_err = result.lin_abs_err ./ params.thresh_x *100;    % relative error (0%-100%)
end
if run_MLE
    result.MLE_t = squeeze(reshape(MLE_t, [step_number, num_start_conditions, num_second_order_weights]));
    result.MLE_s = squeeze(reshape(MLE_s, [step_number, num_start_conditions, num_second_order_weights]));
    result.MLE_abs_err = result.MLE_t - params.thresh_x;                % absolute error (0-1)
    result.MLE_rel_err = result.MLE_abs_err ./ params.thresh_x *100;    % relative error (0%-100%)
end
if version == 13
    result.stochastic_mode = squeeze(reshape(stochastic_mode, [step_number, num_start_conditions, num_second_order_weights]));
end
end
