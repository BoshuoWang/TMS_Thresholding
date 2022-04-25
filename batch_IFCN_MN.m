function batch_IFCN_MN(Batch_id, num_subj_per_batch, step_number)
T_batch_start = tic;

D = dir('DataStorage/Subjects_*.mat');
[~,ind] = max([D.bytes]);
load( fullfile('DataStorage', D(ind).name), 'Subjects', 'Total_Subject_Count', 'y_thresh');

addpath('Statistical-MEP-Model');   % IO model path
addpath('Functions');               % Function and LUT path

T_subj_end = NaN(1, num_subj_per_batch);

txt_fprintf = 'fprintf(''%sThreshold %7.3f%% MSO, relative difference %7.3f%%. ';
 
if Batch_id*num_subj_per_batch > Total_Subject_Count || Batch_id < 1
    fprintf('Batch number %d outside of range for given batch size of %d and total number of subjects (%d).\n', Batch_id, num_subj_per_batch, Total_Subject_Count);
    return
else
    fprintf('Running batch No. %d, subjects No. %d to %d. \n', Batch_id, (Batch_id -1) * num_subj_per_batch + 1, num_subj_per_batch * Batch_id)
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
    fprintf('\n========================================================================\n')
    fprintf('\tSubject %03d -- %s\n', subj_id, datestr(now));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Subject parameters and threshold
    %   Generate paramameters for subject
    params.subj_parameters = Subjects(subj_id).subj_parameters;
    params.start_amplitude = Subjects(subj_id).start_amplitude;
    params.thresh_x = Subjects(subj_id).thresh_x;
    
   %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %   IFCN 
   %   Version 1: 5/10, stepping size 0.01
   %   Version 2: 5/10, stepping size 0.02
   %   Version 3: 10/20, stepping size 0.01
   %   Version 4: 10/20, stepping size 0.02
   for version = 1:4
       T_algo_start = tic;
       eval(sprintf('IFCN%d(subj_cnt) = IFCNmethod(params, version);', version));
       T_algo_end = toc(T_algo_start);
       fprintf('\t\tIFCN version %d: %02d:%02d.%03d.\n\t\t\t', version, fix(T_algo_end/60) , fix(mod(T_algo_end,60)), fix(mod(T_algo_end*1000,1000)) );
       
       eval(sprintf('%s\\n'', %s, %s, %s );', txt_fprintf, '''''', ...
           sprintf('IFCN%d(subj_cnt).th_est', version), ...
           sprintf('IFCN%d(subj_cnt).rel_err', version)));
   end
   
   %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %   Mills and Nithi, 1996, upper and lower thresholds
   %   Version 1: original, threshold search starts from 20% MSO, 1% MSO step size
   %   Version 2: original, threshold search starts from 20% MSO, 2% MSO step size
   %   Version 3: modified, threhsold search starts from suprathrehsold
   %              amplitude obtained during hotspot search, 1% MSO step size   
   %   Version 4: modified, threhsold search starts from suprathrehsold
   %              amplitude obtained during hotspot search, 2% MSO step size   
   for version = 1:4
       T_algo_start = tic;
       eval(sprintf('MillsNithi%d(subj_cnt) = MillsNithiMethod(params, version);', version));
       T_algo_end = toc(T_algo_start);
       fprintf('\t\tMills-Nithi version %d: %02d:%02d.%03d.\n\t\t\t', version, fix(T_algo_end/60) , fix(mod(T_algo_end,60)), fix(mod(T_algo_end*1000,1000)) );
       
       eval(sprintf('%s\\n'', %s, %s, %s );', txt_fprintf, '''''', ...
           sprintf('MillsNithi%d(subj_cnt).th_est', version), ...
           sprintf('MillsNithi%d(subj_cnt).rel_err', version)));
   end
   
    %% %%%%%%%%%%%%%%%%%%%%%%%
    T_subj_end(subj_cnt) = toc(T_subj_start);
    T_run_min = nansum(T_subj_end)/60;
    T_rem_min_est = num_subj_per_batch * nanmean(T_subj_end)/60 - T_run_min;
    fprintf('\tSubject %03d:%3d:%02d:%06.3f.\tRun time: %3d hours %04.1f minutes.\n', subj_id, ...
            floor(T_subj_end(subj_cnt)/3600), floor(mod(T_subj_end(subj_cnt)/60, 60)), mod(T_subj_end(subj_cnt),60), ...
            floor(T_run_min/60), mod(T_run_min,60));
    fprintf('\tAverage: %3d:%02d:%06.3f.\n',...
            floor(nanmean(T_subj_end)/3600),  floor(mod(nanmean(T_subj_end)/60, 60)),  mod(nanmean(T_subj_end),60));
    if subj_cnt > 1
        fprintf('\tRemaining subjects: %3d.\tRemaining time: %3d hours %04.1f minutes (estimated).\n', subj_cnt-1, ...
            floor(T_rem_min_est/60), mod(T_rem_min_est,60) );
    end
end
    
T_batch_end = toc(T_batch_start);
fprintf('\n\nTotal computation time: %2dD %02d:%02d:%06.3f.\n\nSaving data...', ...
        floor(T_batch_end/3600/24), floor(mod(T_batch_end/3600,24)), floor(mod(T_batch_end/60,60)), mod(T_batch_end,60));
    
data_file_name = sprintf('IFCN_MN_batch_No%d.mat', Batch_id);
folder_name = fullfile('DataStorage',sprintf('IFCN_MN_%dsubjperbatch', num_subj_per_batch) );
if ~exist(folder_name,'dir')
    mkdir(folder_name);           % Data storage path
end
save( fullfile(folder_name, data_file_name), 'IFCN*', 'MillsNithi*' );
fprintf('Saved.\n');

end
