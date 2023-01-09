%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   IFCN
%   Version 1: 5/10, stepping size 0.01
%   Version 2: 5/10, stepping size 0.02
%   Version 3: 10/20, stepping size 0.01
%   Version 4: 10/20, stepping size 0.02
function result = IFCNmethod(params, version)
if nargin == 1
    version = 1;
end
if mod(version, 2) == 1 
    ifcn_step_size = 0.01;         	% 1% change in stim amp
else
    ifcn_step_size = 0.02;         	% 2% change in stim amp
end
if version <= 2
    ifcn_maxinum = 10;            	% 5/10
else
    ifcn_maxinum = 20;            	% 10/20
end

% preallocation, maximum steps same as other methods (could be exceeded
% during thresholding
amplitude_list = NaN(1, params.step_number);   
response_bin = false(1, params.step_number); 

%   Integrate Hotspot Info (because hotspot surely above)
step_cnt = 1;                                                               % First step is dummy --> actual search starts with 2 (because hotspot surely above)
amplitude_list(step_cnt) = min( round(params.start_amplitude, 2), 1.3 );    % IFCN method uses integers for amplitudes (1% MSO) precision; maximum amplitude 130% MSO
amplitude_list(step_cnt) = amplitude_list(step_cnt) + ifcn_step_size;       % Dummy step is one step size above, will be reduced to actual start amplitude right as method starts
response_bin(step_cnt) = true;
IFCN_done = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Main
T_start = tic;
while not(IFCN_done)
    % Adjusting the amplitude (e.g., -2% MSO), but should not go to zero
    amplitude_list(step_cnt+1) = max(amplitude_list(step_cnt) - ifcn_step_size, ifcn_step_size);
    
    IFCN_local_positives = 0;
    local_step_cnt = 0;         % Counts at current amplitude
    while (local_step_cnt < ifcn_maxinum)       % maximum 10 stimulation
        step_cnt = step_cnt + 1;
        local_step_cnt = local_step_cnt + 1;
        
        response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt), params.subj_parameters) >= params.y_thresh);   % Directly to binary
        IFCN_local_positives = IFCN_local_positives + response_bin(step_cnt);
        
        if (IFCN_local_positives >= 0.5 * ifcn_maxinum)                         % suprathreshold, e.g., reached 5 out of 10
            break
        elseif ( (local_step_cnt - IFCN_local_positives) > 0.5 * ifcn_maxinum)  % certainly sub-threshold, e.g, already more than 5 of 10 below threshold
            IFCN_done = true;           % Termination
            result.run_time = toc(T_start);
            result.th_est = amplitude_list(step_cnt - local_step_cnt);      % last x(<10) did not lead to 5/10
            result.abs_err = result.th_est - params.thresh_x;               % absolute error (0-1)
            result.rel_err = result.abs_err ./ params.thresh_x *100;
            result.steps = step_cnt - 1;
            result.amplitude_list = amplitude_list(2 : step_cnt);           % Cut it to correct size
            result.response_bin = response_bin(2 : step_cnt);
            break
        else
            amplitude_list(step_cnt+1) = amplitude_list(step_cnt);          % otherwise repeat stimulation
        end
        
    end
end
end
