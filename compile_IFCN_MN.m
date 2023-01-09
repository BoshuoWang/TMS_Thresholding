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
file_name = fullfile(data_folder,sprintf('Compiled_%dsubjs_IFCN_MN.mat', num_subj_per_batch*num_files));
save(file_name, 'IFCN*', 'MillsNithi*');
fprintf('Saved data to: %s\n', file_name);
end