function compile_Bayesian(num_subj_per_batch, num_files)
data_folder = 'DataStorage';

fprintf('Retrieving Bayesian data from: %s\n',data_folder);

if nargin == 0
    num_subj_per_batch = 100;
    num_files = 250;
end
step_number = 200;

folder_name = fullfile(data_folder,sprintf('Bayesian_%dsubjperbatch', num_subj_per_batch) );

for Batch_id = num_files : -1 : 1
    ind = (1 : num_subj_per_batch) + (Batch_id -1) * num_subj_per_batch ;

    for version = 1:3
        data_file_name = sprintf('Bayesian_batch_v%d_No%d.mat', version, Batch_id);
        Data = load( fullfile(folder_name, data_file_name), 'kontsevich*');

        eval(sprintf('kontsevich%d_th_est(ind, 1 : step_number) = cell2mat({Data.kontsevich%d.a_estim}'');', version, version));
        eval(sprintf('kontsevich%d_th_abs_err(ind, 1 : step_number) = cell2mat({Data.kontsevich%d.abs_err}'');', version, version));
        eval(sprintf('kontsevich%d_th_rel_err(ind, 1 : step_number) = cell2mat({Data.kontsevich%d.rel_err}'');', version, version));
        eval(sprintf('kontsevich%d_run_time(ind, 1 : step_number) = cell2mat({Data.kontsevich%d.run_time}'');', version, version));
    end
end
%%
file_name = fullfile(data_folder,sprintf('Compiled_%dsubjs_Bayesian.mat', num_subj_per_batch*num_files));
if exist(file_name, 'file')
    save(file_name,  'kontsevich*', '-append');
else
    save(file_name,  'kontsevich*');
end
fprintf('Saved data to: %s\n',file_name);
end
