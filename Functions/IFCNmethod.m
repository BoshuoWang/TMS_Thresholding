function result = IFCNmethod(params, version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   IFCN
%   Version 1:  5/10, stepping size 0.01 (1% MSO)
%   Version 2:  5/10, stepping size 0.02 (2% MSO)
%   Version 3: 10/20, stepping size 0.01 (1% MSO)
%   Version 4: 10/20, stepping size 0.02 (2% MSO)
if nargin == 1
    version = 1;
end
if mod(version, 2) == 1 
    ifcn_step_size = 0.01;         	% 1% MSO change in stim amp
else
    ifcn_step_size = 0.02;         	% 2% MSO change in stim amp
end
if version <= 2
    ifcn_maxinum = 10;            	%  5/10
else
    ifcn_maxinum = 20;            	% 10/20
end

% preallocation, maximum steps same as other methods
amplitude_list = Inf(1, params.step_number);   
response_bin = false(1, params.step_number); 

%   Integrate hotspot info
step_cnt = 1;
amplitude_list(step_cnt) = round(params.start_amplitude, 2)+ ifcn_step_size;	% IFCN method can only get to 1% precision
response_bin(step_cnt) = true;  % hotspot surely above suprathreshold
IFCN_done = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Main
while not(IFCN_done)
    amplitude_list(step_cnt+1) = amplitude_list(step_cnt) - ifcn_step_size;     % adjusting the amplitude (e.g., -2%)
    
    IFCN_local_positives = 0;
    local_step_cnt = 0;                         % counts at current amplitude
    
    while (local_step_cnt < ifcn_maxinum)       % maximum ifcn_maxinum stimulation
        step_cnt = step_cnt + 1;
        local_step_cnt = local_step_cnt + 1;
        
        % Directly to binary: 0, subthreshold; 1, suprathreshold
        response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt), params.subj_parameters) >= params.y_thresh);
        IFCN_local_positives = IFCN_local_positives + response_bin(step_cnt);
        
        if (IFCN_local_positives >= 0.5 * ifcn_maxinum)                         % suprathreshold, e.g., reached 5 out of 10
            break
        elseif ( (local_step_cnt - IFCN_local_positives) > 0.5* ifcn_maxinum)   % certainly sub-threshold, e.g.,already more than 5 of 10 below threshold
            IFCN_done = true;           % Termination
            result.th_est = amplitude_list(step_cnt - local_step_cnt);          % last x
            result.abs_err = result.th_est - params.thresh_x;                   % absolute error (0-1)
            result.rel_err = result.abs_err ./ params.thresh_x *100;            % relative error (0%-100%)
            result.steps = step_cnt - 1;
            result.amplitude_list = amplitude_list(2 : step_cnt);               % cut to correct size
            result.response_bin = response_bin(2 : step_cnt);
            break
        else
            amplitude_list(step_cnt+1) = amplitude_list(step_cnt);              % otherwise repeat stimulation
        end
    end
end

end
