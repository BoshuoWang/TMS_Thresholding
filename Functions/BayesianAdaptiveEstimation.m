function result = BayesianAdaptiveEstimation(params, version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kontsevich and Tyler, 1999
%   Version 1: Original
%   Version 2: Remove prior information compared original, i.e., uniform pt_ab
%   Version 3: Use estimate as next step.

load('LUT.mat',...
     'size_x', 'size_a', 'size_b',...
     'x_vec',  'a_vec',  'b_vec',...
     'psi_tensor', 'pt_ab_init', 'pt_ab_uniform');

step_number = params.step_number;
amplitude_list = NaN(1, step_number);
response_list = NaN(1, step_number);
run_time      = NaN(1, step_number);

a_estim = NaN(1, step_number);
a_var = NaN(1, step_number);
b_estim = NaN(1, step_number);

amplitude_list(1) = min(params.start_amplitude, max(x_vec));
if version ~= 2
    pt_ab = pt_ab_init;     % Prior information as initial value.
else
    pt_ab = pt_ab_uniform;  % Version 2: uniform pt_ab
end

for step_cnt = 1 : step_number
    T_start = tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Step 1: Calculate probability of getting response r after
    %   presenting test x at the next trial.
    pt_ab_repx = permute(repmat(pt_ab, [1, 1, size_x] ), [3,1,2]);   % building p((a,b)) matrix for all x; faster calculation without loops
    % [x, a, b] <-  [a, b, x] <- [a, b]
    % psi_tensor: p(r|x,θ)
    % pt_ab: p_i (θ)
    % pt_r_SUB_x: ∫_(R_+^n) [p_i (θ)p(r|x,θ) d^n θ]
    pt_r_SUB_x(:, 2) = sum( psi_tensor .* pt_ab_repx, [2,3] );	% Suprathreshold response, avoid preallocation by assigning second row first
    % p(suprathreshold | x) = sum_(a,b) p(suprathreshold,(a,b) | x) = sum(a,b) ( p(suprathreshold | x, (a,b)) * p((a,b)) )
    %  [x, 2]                      [x, a, b] .* [x, a, b]
    pt_r_SUB_x(:, 1) = 1 - pt_r_SUB_x(:, 2);            % Subthreshold response
    % p(subthreshold | x) = 1 - p(suprathreshold | x)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Step 2: Estimate by Bayes' rule the posterior probability of each
    %   psychometric function given that the next trial will produce the
    %   response r to the test of the intensity x.
    pt_ab_SUB_xr(:, :, :, 2) = permute(       psi_tensor  .* pt_ab_repx ./ repmat(pt_r_SUB_x(:, 2), [1, size_a, size_b]) , [2,3,1]);
    % pt_ab_SUB_xr: p_i (θ|x,r)
    % p((a,b) | suprathreshold, x) =      p(suprathreshold | x, (a,b)) *    p((a,b))       /     p_t(suprathreshold | x))
    % [a, b, x, 2]                  [a, b, x] <---       [x, a, b]    .* [x, a, b]        ./ [x, a, b]    <- [x,1]
    pt_ab_SUB_xr(:, :, :, 1) = permute(  (1 - psi_tensor) .* pt_ab_repx ./ repmat(pt_r_SUB_x(:, 1), [1, size_a, size_b]) , [2,3,1]);
    pt_ab_SUB_xr = max(pt_ab_SUB_xr, eps);   %   Intermediate Step (added by S.goetz 06/14/2010, shall prevent ld0); mod. by BWang, 3/10/2021
        
    if (version ~= 3)   % Original stepping methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Step 3: Estimate the entropy of the probability density
        %   function over the space of psychometric functions, given
        %   that at the next trial a test of intensity x will produce
        %   the response r.    Conditional entropy H(a,b | r,x)
        H_xr = -1 * squeeze(sum(pt_ab_SUB_xr .* log2(pt_ab_SUB_xr), [1, 2]));
        % [x, 2]  <-----         [a, b, x, 2]
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Step 4: Estimate the expected entropy for each test intensity x.
        EH_x =  sum(H_xr .* pt_r_SUB_x, 2); 	% Expected entropy
        % [x]  ;  [x, 2] .* [x, 2]
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Step 5: Find the test intensity that has the minimum expected
        %   entropy.
        [~, x_index] = min(EH_x);          % EH_min_index means index of stimulation amplitude
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Alternative: use threshold estimate from previous step for
        %   current step. For first step, use starting amplitude
        if (step_cnt > 1)
            [~, x_index] = min( abs(a_estim(step_cnt-1) - x_vec) );     % Find corresponding amplitude in look-up table
        else
            [~, x_index] = min( abs(amplitude_list(1) - x_vec) );       % Find amplitude of first trial in look-up table
        end
    end
    amplitude_list(step_cnt) = x_vec(x_index);	% Transfer index to stimulation amplitude value
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Step 6: Run a trial with intensity x_(t+1) to obtain the response r_(t+1).
    %   Stimulate either with start amplitude (t = 0, first step) or next (chosen) amplitude (later steps)
    response_list(step_cnt) = virtstimulate( amplitude_list(step_cnt), params.subj_parameters );
    response_bin = 1 + ( response_list(step_cnt) >= params.y_thresh );
    % 0/1 logic to 1/2-indices: 2 = suprathreshold, 1 = subthreshold
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %	Estimate new parameters:
    %   Step 7: Keep the posterior probability distribution from Step 2
    %   that corresponds to the completed trial. 
    pt_ab = pt_ab_SUB_xr(:, :, x_index, response_bin);
    % [a, b]     ; [a, b, x, 2]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Step 8: Find a new estimate of the psychometric function based on
    %   the new posterior probability distribution p_(t+1)_(a,b).
    a_estim(step_cnt) = sum( a_vec * pt_ab  );    % [1, a]*[a, b]
    b_estim(step_cnt) = sum( pt_ab * b_vec' );   % [a, b]*[b, 1]
    %   Error/Deviation, Termination Information
    a_var(step_cnt) = sum( a_vec.^2 * pt_ab ) - a_estim(step_cnt)^2;
    %   Variance from displacement theorem (E{(X-E{X})^2} = E{X^2}-(E{X})^2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Step 9: Return to Step 1, unless a specified number of trials is completed. 
    run_time(step_cnt) = toc(T_start);
   
end
abs_err = a_estim - params.thresh_x; % absolute error (0-1)
rel_err = abs_err ./ params.thresh_x *100;

result = struct('amplitude_list', amplitude_list, ...
                'response_list', response_list, ...
                'a_estim', a_estim, ...
                'a_var', a_var, ...
                'b_estim', b_estim,...
                'abs_err', abs_err,...
                'rel_err', rel_err, ...
                'run_time', run_time ...
               );
end
