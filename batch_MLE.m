function batch_MLE(Batch_id, num_subj_per_batch, step_number, version, fid)
if nargin <= 4
    fid = 1;
end
T_batch_start = tic;

data_folder = 'DataStorage';
if ~exist('Total_Subject_Count', 'var')
    Total_Subject_Count = 25e3;    % 25k
end
file_name = fullfile(data_folder, sprintf('Subjects_%d_reprocessed.mat', Total_Subject_Count));
load( file_name, 'Subjects', 'Total_Subject_Count', 'y_thresh');

addpath('Statistical-MEP-Model');   % IO model path
addpath('Functions');               % Function and LUT path

T_subj_end = NaN(1, num_subj_per_batch);
fprintf(fid, 'Number of steps: %d\n', step_number);
txt_fprintf = 'fprintf(fid, ''%sThreshold %7.3f%% MSO, relative difference %7.3f%%. ';

if Batch_id*num_subj_per_batch > Total_Subject_Count || Batch_id < 1
    fprintf(fid, 'Batch number %d outside of range for given batch size of %d and total number of subjects (%d).\n', Batch_id, num_subj_per_batch, Total_Subject_Count);
    return
else
    fprintf(fid, 'Running batch No. %d, subjects No. %d to %d. \n', Batch_id, (Batch_id -1) * num_subj_per_batch + 1, num_subj_per_batch * Batch_id);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters:
params.step_number = step_number;   % Package parameters for passing into functions
params.y_thresh = y_thresh;
params.opts = optimset('Display', 'off', 'MaxFunEvals', 500000, 'FunValCheck', 'on', 'MaxIter', 10000, 'TolFun', 1e-6, 'TolX', 1e-8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Main Loop:
for subj_cnt = num_subj_per_batch : -1 : 1        % Reverse loop, elimitnate need to pre-allocate space for many variables
    T_subj_start = tic;
    subj_id = subj_cnt + (Batch_id -1) * num_subj_per_batch;
    fprintf(fid, '\n========================================================================\n');
    fprintf(fid, '\tSubject %03d -- %s\n', subj_id, datetime('now'));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Subject parameters and threshold
    %   Generate paramameters for subject
    params.subj_parameters = Subjects(subj_id).subj_parameters;
    params.start_amplitude = Subjects(subj_id).p95_start.amplitude;
    params.thresh_x = Subjects(subj_id).relative_frequency.p50_lin;

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Maximum Likelihood Estimation methods
    %   Different versions, rewritten Nov. 2022
    %   Version 1: MLE-1 estimation & stepping, explicit evaluation (Awiszus, 2003)
    %   Version 2: MLE-1 estimation & stepping, direct maximization
    %   Version 3: MAP-1 estimation & stepping, direct maximization
    %
    %   Version 4: MLE-2 estimation & stepping, explicit evaluation
    %	Version 5: MLE-2 estimation & stepping, direct maximization
    %   Version 6: MAP-2 estimation & stepping, direct maximization
    %
    %	Version 7: MAP-2 estimation, min variance stepping, direct maximization
    %   Version 8: MAP-2 estimation, min inertia stepping, direct maximization

    for MLE_ver = version
        T_algo_start = tic;
        eval(sprintf('MLE%d(subj_cnt) = MLEstimation(params, MLE_ver);', MLE_ver));
        T_algo_end = toc(T_algo_start);
        fprintf(fid, '\t\tMLE version %d: %02d:%02d.%03d.\n\t\t\t', MLE_ver, fix(T_algo_end/60) , fix(mod(T_algo_end,60)), fix(mod(T_algo_end*1000,1000)) );

        eval(sprintf('%s\\n'', %s, %s * 100, %s );', txt_fprintf, '''''', ...
            sprintf('MLE%d(subj_cnt).thresh_est(end)', MLE_ver), ...
            sprintf('MLE%d(subj_cnt).rel_err(end)', MLE_ver)));
    end
          
    %% %%%%%%%%%%%%%%%%%%%%%%%
    T_subj_end(subj_cnt) = toc(T_subj_start);
    T_run_min = sum(T_subj_end, 'omitnan')/60;
    T_subj_end_ave = mean(T_subj_end, 'omitnan');
    T_rem_min_est = num_subj_per_batch * T_subj_end_ave/60 - T_run_min;
    fprintf(fid, '\tSubject %03d:%3d:%02d:%06.3f.\tRun time: %3d hours %04.1f minutes.\n', subj_id, ...
            floor(T_subj_end(subj_cnt)/3600), floor(mod(T_subj_end(subj_cnt)/60, 60)), mod(T_subj_end(subj_cnt), 60), ...
            floor(T_run_min/60), mod(T_run_min,60));
    fprintf(fid, '\tAverage: %3d:%02d:%06.3f.\n',...
            floor(T_subj_end_ave/3600),  floor(mod(T_subj_end_ave/60, 60)),  mod(T_subj_end_ave, 60));
    if subj_cnt > 1
        fprintf(fid, '\tRemaining subjects: %3d.\tRemaining time: %3d hours %04.1f minutes (estimated).\n', subj_cnt-1, ...
            floor(T_rem_min_est/60), mod(T_rem_min_est, 60) );
    end
end
    
T_batch_end = toc(T_batch_start);
fprintf(fid, '\n\nTotal computation time: %2dD %02d:%02d:%06.3f.\n\nSaving data...', ...
        floor(T_batch_end/3600/24), floor(mod(T_batch_end/3600,24)), floor(mod(T_batch_end/60,60)), mod(T_batch_end,60));
    
data_file_name = sprintf('MLE_v%d_No%d.mat', version, Batch_id);
folder_name = fullfile('DataStorage',sprintf('MLE_%dsubjperbatch', num_subj_per_batch) );
if ~exist(folder_name,'dir')
    mkdir(folder_name);           % Data storage path
end

save( fullfile(folder_name, data_file_name), 'MLE*');
fprintf(fid, 'Saved.\n');

end
