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
%   the range of [1, 130]. 
%   Need to convert to [0, 1] before passing to virtstimulate function and
%   returning results. 

if nargin == 1
    version = 1;
end

switch version
    case {1,3}
        NM_step_size = 1;   % 1% MSO change in stim amp 
    case {2,4}
        NM_step_size = 2;   % 2% MSO change in stim amp  
end
NM_maximum = 10;         	% 10/10

% preallocation, maximum steps same as other methods (could be exceeded
% during thresholding
amplitude_list = NaN(1, params.step_number);   
response_bin = false(1, params.step_number);
response_counts = zeros(130, 2); % Negative counts and positive counts, for 1%-130% MSO

%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Main
T_start = tic;
step_cnt = 1;
if version == 1 || version == 2              % original method, threshold search starts from 20% MSO
    amplitude_list(step_cnt) = 20;           % First stim at 20%
    response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);       % Directly to binary
    response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;        % convert 0/1 response logic to 1/2 indices for neg/pos counts

    % Find first suptrathreshold response
    while ~response_bin(step_cnt)           
        step_cnt = step_cnt + 1;
        amplitude_list(step_cnt) = min(amplitude_list(step_cnt-1) + 10, 130);   % Increasing by 10% until response found, maximum amplitude 130% MSO
        response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);   % Directly to binary
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
            response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;    % convert 0/1 response logic to 1/2 indices for neg/pos counts
    end

    % Find lower threshold
    amplitude_list(step_cnt+1) = max(amplitude_list(step_cnt) - NM_step_size, NM_step_size);
    negative_count = response_counts(amplitude_list(step_cnt+1), 1);
    while (negative_count < NM_maximum)       % maximum 10 stimulation with negative response
        step_cnt = step_cnt + 1;
        response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);   % Directly to binary
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
            response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;    % convert 0/1 response logic to 1/2 indices for neg/pos counts

        if  response_bin(step_cnt)  % positive response, decrease amplitude
            amplitude_list(step_cnt+1) = max(amplitude_list(step_cnt) - NM_step_size, NM_step_size);
        else                        % negative response, repeat until either have 10 negative responses or 1 positive response
            amplitude_list(step_cnt+1) = amplitude_list(step_cnt);
        end
        negative_count = response_counts(amplitude_list(step_cnt+1), 1);
    end
    result.lo_th = amplitude_list(step_cnt)/100;

    % Find upper threshold
    amplitude_list(step_cnt+1) = find(response_counts(:, 1)~=0, 1, 'last') + 1; % start above largest amplitude that had negative response
    positive_count = response_counts(amplitude_list(step_cnt+1), 2);
    while (positive_count  < NM_maximum)       % maximum 10 stimulation with positive response
        step_cnt = step_cnt + 1;
        response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);   % Directly to binary
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
            response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;    % convert 0/1 response logic to 1/2 indices for neg/pos counts

        if  response_bin(step_cnt)  % positive response, repeat until either have 10 positive responses or 1 negative response
            amplitude_list(step_cnt+1) = amplitude_list(step_cnt);
        else                        % negative response, increase amplitude
            amplitude_list(step_cnt+1) = min(amplitude_list(step_cnt) +  NM_step_size, 130);
        end
        positive_count = response_counts(amplitude_list(step_cnt+1), 2);
    end
    result.up_th = amplitude_list(step_cnt)/100;
else
    amplitude_list(step_cnt) = min(round(params.start_amplitude * 100), 130);       % can only get to 1% precision, maximum amplitude 130% MSO
    response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);       % Directly to binary
    response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;        % convert 0/1 response logic to 1/2 indices for neg/pos counts

    % Find first subthreshold response
    while  response_bin(step_cnt)           
        step_cnt = step_cnt + 1;
        amplitude_list(step_cnt) = max(amplitude_list(step_cnt-1) - 10, NM_step_size);   % Decreasing by 10% until no response found, minimum amplitude should not reach zeros
        response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);   % Directly to binary
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
            response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;    % convert 0/1 response logic to 1/2 indices for neg/pos counts
    end

    % Find upper threshold
    amplitude_list(step_cnt+1) = min(amplitude_list(step_cnt) + NM_step_size, 130);
    positive_count = response_counts(amplitude_list(step_cnt+1), 2);
    while (positive_count < NM_maximum)       % maximum 10 stimulation with positive response
        step_cnt = step_cnt + 1;
        response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);   % Directly to binary
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
            response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;    % convert 0/1 response logic to 1/2 indices for neg/pos counts

        if ~response_bin(step_cnt)  % negative response, increase amplitude
            amplitude_list(step_cnt+1) = min(amplitude_list(step_cnt) + NM_step_size, 130);
        else                        % positive response, repeat until either have 10 positive responses or 1 negative response
            amplitude_list(step_cnt+1) = amplitude_list(step_cnt);
        end
        positive_count = response_counts(amplitude_list(step_cnt+1), 2);
    end
    result.up_th = amplitude_list(step_cnt)/100;

    % Find lower threshold
    amplitude_list(step_cnt+1) = find(response_counts(:, 2)~=0, 1, 'first') - 1; % start below smallest amplitude that had positive response
    negative_count = response_counts(amplitude_list(step_cnt+1), 1);
    while (negative_count  < NM_maximum)       % maximum 10 stimulation with negative response
        step_cnt = step_cnt + 1;
        response_bin(step_cnt) = logical(virtstimulate(amplitude_list(step_cnt)/100, params.subj_parameters) >= params.y_thresh);   % Directly to binary
        response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) = ...
            response_counts(amplitude_list(step_cnt), response_bin(step_cnt)+1) + 1;    % convert 0/1 response logic to 1/2 indices for neg/pos counts

        if ~response_bin(step_cnt)  % negative response, repeat until either have 10 negative responses or 1 positive response
            amplitude_list(step_cnt+1) = amplitude_list(step_cnt);
        else                        % positive response, decrease amplitude
            amplitude_list(step_cnt+1) = max(amplitude_list(step_cnt) -  NM_step_size, NM_step_size);
        end
        negative_count = response_counts(amplitude_list(step_cnt+1), 1);
    end
    result.lo_th = amplitude_list(step_cnt)/100;
end

result.run_time = toc(T_start);
result.th_est = (result.up_th + result.lo_th)/2;
result.abs_err = result.th_est - params.thresh_x; % absolute error (0-1)
result.rel_err = result.abs_err ./ params.thresh_x *100;
result.steps = step_cnt;
result.amplitude_list = amplitude_list(1 : step_cnt)/100;    % Cut it to correct size and scale to [0,1]
result.response_bin = response_bin(1 : step_cnt);
            
end