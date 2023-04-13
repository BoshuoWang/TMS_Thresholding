function compile_IFCN_MN(num_subj_per_batch, num_files)
data_folder = 'DataStorage';

fprintf('Retrieving IFCN and MN data from: %s\n', data_folder);

if nargin == 0
    num_subj_per_batch = 100;
    num_files = 250;
end

folder_name = fullfile('DataStorage',sprintf('IFCN_MN_%dsubjperbatch', num_subj_per_batch) );

for Batch_id = 1 : num_files
    ind = (1 : num_subj_per_batch) + (Batch_id - 1) * num_subj_per_batch ;

    data_file_name = sprintf('IFCN_MN_batch_No%d.mat', Batch_id);
    Data = load( fullfile(folder_name, data_file_name), 'IFCN*', 'MillsNithi*' );
    for jj = 1:4
        eval(sprintf('IFCN%d.steps(ind, 1)  =  	[Data.IFCN%d.steps]'';', jj, jj));
        eval(sprintf('IFCN%d.th_est(ind, 1) =  	[Data.IFCN%d.th_est]'';', jj, jj));
        eval(sprintf('IFCN%d.abs_err(ind, 1) =  [Data.IFCN%d.abs_err]'';', jj, jj));
        eval(sprintf('IFCN%d.rel_err(ind, 1) =  [Data.IFCN%d.rel_err]'';', jj, jj));
        eval(sprintf('IFCN%d.run_time(ind, 1) = [Data.IFCN%d.run_time]'';', jj, jj));

        eval(sprintf('MillsNithi%d.steps(ind, 1)    = [Data.MillsNithi%d.steps]'';', jj, jj));
        eval(sprintf('MillsNithi%d.th_est(ind, 1)   = [Data.MillsNithi%d.th_est]'';', jj, jj));
        eval(sprintf('MillsNithi%d.abs_err(ind, 1)  = [Data.MillsNithi%d.abs_err]'';', jj, jj));
        eval(sprintf('MillsNithi%d.rel_err(ind, 1)  = [Data.MillsNithi%d.rel_err]'';', jj, jj));
        eval(sprintf('MillsNithi%d.run_time(ind, 1) = [Data.MillsNithi%d.run_time]'';', jj, jj));
    end
end
%%
data_type = {'IFCN3', 'IFCN4', 'IFCN1', 'IFCN2', 'MillsNithi1', 'MillsNithi2', 'MillsNithi3', 'MillsNithi4'};
err_str = {'Absolute', 'Relative'} ;
for ll = 1 : length(data_type)
    fprintf('Method: %s. \n', data_type{ll});
    procedure_steps = Data.(data_type{ll}).steps;
    fprintf('Mean and std of steps: %3.1f +/- %3.1f. Median steps: %3.1f\n', mean(procedure_steps), std(procedure_steps), median(procedure_steps));
        
    for kk = 1 : length(err_str) 
        fprintf('Error type: %s. \n', err_str{kk});
        if kk == 1
            th_err = Data.(data_type{ll}).(sprintf('%s_err', err_type))*100;
        else
            th_err = Data.(data_type{ll}).(sprintf('%s_err', err_type));
        end
        stats = regstats(th_err, procedure_steps, 'linear');
        Data.(data_type{ll}).(sprintf('%s_err_stats', err_type)) = stats;
        fprintf('Mean and std of error: %3.2f%% +/- %3.2f%%. Median error: %3.2f%%\n', mean(th_err),std(th_err),median(th_err));
        fprintf('Slope: %3.3f%% .\n\n', stats.beta(2));        
    end
end

file_name = fullfile(data_folder,sprintf('Compiled_%dsubjs_IFCN_MN.mat', num_subj_per_batch*num_files));
save(file_name, 'IFCN*', 'MillsNithi*');
fprintf('Saved data to: %s\n', file_name);
end