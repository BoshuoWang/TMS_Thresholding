function result = MillsNithiMethod(params, version)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Mills and Nithi, 1996, upper and lower thresholds
%   Version 1: original, threshold search starts from 20% MSO, 1% MSO step size
%   Version 2: original, threshold search starts from 20% MSO, 2% MSO step size
%   Version 3: modified, threhsold search starts from suprathrehsold
%              amplitude obtained during hotspot search, 1% MSO step size
%   Version 4: modified, threhsold search starts from suprathrehsold
%              amplitude obtained during hotspot search, 2% MSO step size

%   For easier indexing, the amplitudes in this function are integers in
%   the range of [1,100]. 
%   Need to convert to [0,1] before passing to virtstimulate function and
%   returning results. 

if nargin == 1
    version = 1;
end

switch version
    case {1,3}
        NM_step_size = 1;   % 1% change in stim amp 
    case {2,4}
        NM_step_size = 2;   % 2% change in stim amp  
end
NM_maximum = 10;         	% 5/10

amplitude_list = Inf(1, params.step_number);   
response_bin = false(1, params.step_number);
response_counts = zeros(200, 2); % negative counts and positive counts, for 1%-200% MSO

%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Main
step_cnt = 1;
if version == 1 || version == 2
    amplitude_list(step_cnt) = 20;           % First stim at 20%
    response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);
    response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;        % convert 0/1 response logic to 1/2 indices for neg/pos counts

    % Find first suptrathreshold response
    while ~response_bin(step_cnt)           
        step_cnt = step_cnt + 1;
        amplitude_list(step_cnt) = amplitude_list(step_cnt-1) + 10;                     % increasing by 10% until response found
        response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
            response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;    % convert 0/1 response logic to 1/2 indices for neg/pos counts
    end

    % Find lower threshold
    amplitude_list(step_cnt+1) = amplitude_list(step_cnt) - NM_step_size;
    positive_count = response_counts(amplitude_list(step_cnt+1), 1);
    while (positive_count < NM_maximum)                                                 % maximum 10 stimulation
        step_cnt = step_cnt + 1;
        response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
            response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;    % convert 0/1 response logic to 1/2 indices for neg/pos counts

        if response_bin(step_cnt)                                                       % positive response
            amplitude_list(step_cnt+1) = amplitude_list(step_cnt) - NM_step_size;
        else 
            amplitude_list(step_cnt+1) = amplitude_list(step_cnt);
        end
        positive_count = response_counts(amplitude_list(step_cnt+1), 1);
    end
    result.lo_th = amplitude_list(step_cnt)/100;

    % Find upper threshold
    amplitude_list(step_cnt+1) = find(response_counts(:, 1)~=0, 1, 'last') + 1;
    negative_count = response_counts(amplitude_list(step_cnt+1), 2);
    while (negative_count  < NM_maximum)       % maximum 10 stimulation
        step_cnt = step_cnt + 1;
        response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
            response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;    % convert 0/1 response logic to 1/2 indices for neg/pos counts

        if response_bin(step_cnt)                                                       % positive response
            amplitude_list(step_cnt+1) = amplitude_list(step_cnt);
        else 
            amplitude_list(step_cnt+1) = amplitude_list(step_cnt) +  NM_step_size;
        end
        negative_count = response_counts(amplitude_list(step_cnt+1), 2);
    end
    result.up_th = amplitude_list(step_cnt)/100;

else
    amplitude_list(step_cnt) = round(params.start_amplitude * 100);                     % can only get to 1% precision
    
    response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);
    response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;        % convert 0/1 response logic to 1/2 indices for neg/pos counts

    % Find first subthreshold response
    while response_bin(step_cnt)           
        step_cnt = step_cnt + 1;
        amplitude_list(step_cnt) = amplitude_list(step_cnt-1) - 10;                     % decreasing by 10% until response found
        response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
            response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;    % convert 0/1 response logic to 1/2 indices for neg/pos counts
    end

    % Find upper threshold
    amplitude_list(step_cnt+1) = amplitude_list(step_cnt) + NM_step_size;
    positive_count = response_counts(amplitude_list(step_cnt+1), 2);
    while (positive_count < NM_maximum)                                                 % maximum 10 stimulation
        step_cnt = step_cnt + 1;
        response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
            response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;    % convert 0/1 response logic to 1/2 indices for neg/pos counts

        if ~response_bin(step_cnt)                                                      % negative response
            amplitude_list(step_cnt+1) = amplitude_list(step_cnt) + NM_step_size;
        else 
            amplitude_list(step_cnt+1) = amplitude_list(step_cnt);
        end
        positive_count = response_counts(amplitude_list(step_cnt+1), 2);
    end
    result.up_th = amplitude_list(step_cnt)/100;

    % Find lower threshold
    amplitude_list(step_cnt+1) = find(response_counts(:, 2)~=0, 1, 'first') - 1;
    negative_count = response_counts(amplitude_list(step_cnt+1), 1);
    while (negative_count  < NM_maximum)                                                % maximum 10 stimulation
        step_cnt = step_cnt + 1;
        response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
            response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;    % convert 0/1 response logic to 1/2 indices for neg/pos counts

        if ~response_bin(step_cnt)                                                      % negative response
            amplitude_list(step_cnt+1) = amplitude_list(step_cnt);
        else 
            amplitude_list(step_cnt+1) = amplitude_list(step_cnt) -  NM_step_size;
        end
        negative_count = response_counts(amplitude_list(step_cnt+1), 1);
    end
    result.lo_th = amplitude_list(step_cnt)/100;
end

result.th_est = (result.up_th + result.lo_th)/2;
result.abs_err = result.th_est - params.thresh_x;           % absolute error (0-1)
result.rel_err = result.abs_err ./ params.thresh_x *100;    % relative error (0%-100%)
result.steps = step_cnt;
result.amplitude_list = amplitude_list(1 : step_cnt)/100;	% cut to correct size
result.response_bin = response_bin(1 : step_cnt);
            
end
