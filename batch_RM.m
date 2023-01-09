function batch_RM(Batch_id, num_subj_per_batch, step_number, version, fid)
if nargin < 4
    fid = 1;
end

data_file_name = sprintf('RobMon_v%d_No%d.mat', version, Batch_id);
folder_name = fullfile('DataStorage',sprintf('RobMon_%dsubjperbatch', num_subj_per_batch) );
if ~exist(folder_name, 'dir')
    mkdir(folder_name);           % Data storage path
end
if exist(fullfile(folder_name, data_file_name), 'file')
    return
end

%%
T_batch_start = tic;

D = dir('DataStorage/Subjects_*.mat');
[~,ind] = max([D.bytes]);
load( fullfile('DataStorage', D(ind).name), 'Subjects', 'Total_Subject_Count', 'y_thresh');

addpath('Statistical-MEP-Model');   % IO model path
addpath('Functions');               % Function and LUT path

T_subj_end = NaN(1, num_subj_per_batch);

fprintf(fid, 'Number of steps: %d\n', step_number);
 
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
    params.start_amplitude = Subjects(subj_id).start_amplitude;
    params.thresh_x = Subjects(subj_id).thresh_x;
        
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Stochastic Approximation, Robbins-Monro, fully log, all versions
    %   Version 1, 2:       non-adaptive 1/n, 1st & 2nd order (original versions 9, 10)
    %   Version 3, 4:       adaptive 1/n, 1st & 2nd order (original versions 1/11/81, 2/32/42)
    %   Version 5, 6:       adaptive 2^(-n), 1st & 2nd order (original versions 3, 4)
    %   Version 7, 8:       adaptive 2^(-n/2), 1st & 2nd order (original versions 5/15, 6)
    %   Version 9, 10:      adaptive 2^(-n/2), 1st & 2nd order, de- and increasing steps (original versions 7, 8)
    %   Version 11, 12:     adaptive 1/n, 1st & 2nd order, step increases if plateau reached (original version 21)
    %   Version 13:         Stochastic Newton, starts as RM (adaptive 1/n, 1st order) until robust estimation with linear regression (orginal version snewta)
    %   Version 7-3:        adaptive 2^(-n/2) switches to 1/n, 1st order i.e., rapid approach switching to a.s. convergence. Transition: max step stize < +/- 0.015 = 1.5% MSO (original version 151)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    is_analog = true;
            
    switch version
        case {1, 5, 7, 9, 73, 11}
            num_start_conditions = 21;
            
            T_algo_start = tic;
            eval(sprintf('RobMonA%d_s(subj_cnt) = StochasticApproximation(params,  version, is_analog, num_start_conditions);', version));
            T_algo_end = toc(T_algo_start);
            fprintf(fid, '\t\tRobbins-Monro analog version %d (parameter sweep): %02d:%02d.%03d.\n', version, fix(T_algo_end/60) , fix(mod(T_algo_end,60)), fix(mod(T_algo_end*1000,1000)) );
            
            if (version < 10)
                is_analog = false;
                
                T_algo_start = tic;
                eval(sprintf('RobMonD%d_s(subj_cnt) = StochasticApproximation(params,  version, is_analog, num_start_conditions);', version));
                T_algo_end = toc(T_algo_start);
                fprintf(fid, '\t\tRobbins-Monro digital version %d (parameter sweep): %02d:%02d.%03d.\n', version, fix(T_algo_end/60) , fix(mod(T_algo_end,60)), fix(mod(T_algo_end*1000,1000)) );
            end
            
        case {3}
            num_start_conditions = 1;
            num_second_order_weights = 1;
            run_lin = true;
            run_MLE = true;
            
            T_algo_start = tic;
            RobMonA3(subj_cnt) = StochasticApproximation(params, version, is_analog, num_start_conditions, num_second_order_weights, run_lin, run_MLE);
            T_algo_end = toc(T_algo_start);
            fprintf(fid, '\t\tRobbins-Monro analog version %d: %02d:%02d.%03d.\n', version, fix(T_algo_end/60) , fix(mod(T_algo_end,60)), fix(mod(T_algo_end*1000,1000)) );
            
            is_analog = false;
            num_start_conditions = 21;
            T_algo_start = tic;
            eval(sprintf('RobMonD%d_s(subj_cnt) = StochasticApproximation(params,  version, is_analog, num_start_conditions);', version));
            T_algo_end = toc(T_algo_start);
            fprintf(fid, '\t\tRobbins-Monro digital version %d (parameter sweep): %02d:%02d.%03d.\n', version, fix(T_algo_end/60) , fix(mod(T_algo_end,60)), fix(mod(T_algo_end*1000,1000)) );
        case {2,6,8}
            T_algo_start = tic;
            eval(sprintf('RobMonA%d(subj_cnt) = StochasticApproximation(params,  version, is_analog);', version));
            T_algo_end = toc(T_algo_start);
            fprintf(fid, '\t\tRobbins-Monro analog version %d: %02d:%02d.%03d.\n', version, fix(T_algo_end/60) , fix(mod(T_algo_end,60)), fix(mod(T_algo_end*1000,1000)) );
           
        case {4}
            % Original version roba42: parameter search for control
            % sequence start point and second order weight
            num_start_conditions = 21;
            num_second_order_weights = 41;
            
            T_algo_start = tic;
            RobMonA4_sw(subj_cnt) = StochasticApproximation(params, version, is_analog, num_start_conditions, num_second_order_weights);
            T_algo_end = toc(T_algo_start);
            fprintf(fid, '\t\tRobbins-Monro analog version %d (2D parameter sweep): %02d:%02d.%03d.\n', version, fix(T_algo_end/60) , fix(mod(T_algo_end,60)), fix(mod(T_algo_end*1000,1000)) );
            
        case {13}
            num_start_conditions = 1;
            num_second_order_weights = 1;
            run_lin = true;
            run_MLE = true;
            
            T_algo_start = tic;
            RobMonA13(subj_cnt) = StochasticApproximation(params, version, is_analog, num_start_conditions, num_second_order_weights, run_lin, run_MLE);
            T_algo_end = toc(T_algo_start);
            fprintf(fid, '\t\tRobbins-Monro analog version %d (stochastic Newton): %02d:%02d.%03d.\n', version, fix(T_algo_end/60) , fix(mod(T_algo_end,60)), fix(mod(T_algo_end*1000,1000)) );
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
    
save( fullfile(folder_name, data_file_name), 'RobMon*', '-v7.3');
fprintf(fid, 'Saved.\n');

end
