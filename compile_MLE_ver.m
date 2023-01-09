function compile_MLE_ver(version, num_subj_per_batch, num_files)
data_folder = 'DataStorage';

fprintf('Retrieving MLE data from: %s\n',data_folder);

if nargin == 1
    num_subj_per_batch = 100;
    num_files = 250;
end
step_number = 200;

folder_name = fullfile(data_folder,sprintf('MLE_%dsubjperbatch', num_subj_per_batch) );

for Batch_id = num_files : -1 : 1
    ind = (1 : num_subj_per_batch) + (Batch_id -1) * num_subj_per_batch ;
    
    
    data_file_name = sprintf('MLE_v%d_No%d.mat', version, Batch_id);
    Data = load( fullfile(folder_name, data_file_name), 'MLE*');

    eval(sprintf('MLE%d_th_est(ind, 1 : step_number) = cell2mat({Data.MLE%d.thresh_est}'');', version, version));
    eval(sprintf('MLE%d_th_abs_err(ind, 1 : step_number) = cell2mat({Data.MLE%d.abs_err}'');', version, version));
    eval(sprintf('MLE%d_th_rel_err(ind, 1 : step_number) = cell2mat({Data.MLE%d.rel_err}'');', version, version));
    eval(sprintf('MLE%d_run_time(ind, 1 : step_number) = cell2mat({Data.MLE%d.run_time}'');', version, version));

end
%%
file_name = fullfile(data_folder,sprintf('Compiled_%dsubjs_MLE_%d.mat', num_subj_per_batch*num_files, version));
if exist(file_name, 'file')
    save(file_name,  'MLE*', '-append');
else
    save(file_name,  'MLE*');
end
fprintf('Saved data to: %s\n', file_name);
end
