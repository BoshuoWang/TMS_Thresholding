function compile_reSubjects(num_subj_per_batch, num_files)
data_folder = 'DataStorage';
fprintf('Retrieving reprocessed subject data from: %s\n',data_folder);
if nargin == 0
    num_subj_per_batch = 100;
    num_files = 250;
end
if ~exist('Total_Subject_Count', 'var')
    Total_Subject_Count = 25e3;    % 25k
end
data_file_name = sprintf('Subjects_%d.mat', Total_Subject_Count);
load( fullfile('DataStorage', data_file_name), 'Total_Subject_Count', 'y_thresh');

folder_name = fullfile('DataStorage', sprintf('reSubjects_%dsubjperbatch', num_subj_per_batch) );
for Batch_id = num_files : -1 : 1
    data_file_name = sprintf('reSubjects_No%d.mat', batch_id);
    Data = load( fullfile(folder_name, data_file_name));   
    
    ind = (1 : num_subj_per_batch) + (Batch_id -1) * num_subj_per_batch ;
    Subjects(ind) = Data.NewSubjects;
end
%%
file_name = fullfile(data_folder, sprintf('Subjects_%d_reprocessed.mat', Total_Subject_Count));
save(file_name, 'Subjects', 'Total_Subject_Count', 'y_thresh');
fprintf('Saved data to: %s\n', file_name);
end
