addpath('Statistical-MEP-Model');   % IO model path

data_folder = 'DataStorage';
if ~exist('Total_Subject_Count', 'var')
    Total_Subject_Count = 25e3;    % 25k
end
file_name = fullfile(data_folder, sprintf('Subjects_%d_reprocessed_v1.mat', Total_Subject_Count));
load(file_name, 'Subjects', 'Total_Subject_Count', 'y_thresh');

%%
stim_amp_eps = 0.001e-2;                            %
stim_amp_step = 0.2e-2;                             % 0.2% MSO
stim_amp_min = stim_amp_step;
stim_amp_max = 1.5;                                 % 150% MSO

noisefree = true;
all_stim_vec = (stim_amp_min : stim_amp_step : stim_amp_max);               % 0%-150% MSO, in steps of 0.2% MSO (750 pulse amplitude)

oldFields = {'old_start', 'NL50uV'};
newFields = {'NL_500uV', 'NL_50uV'};
for ii = 1 : length(oldFields)
    [Subjects.(newFields{ii})] = Subjects.(oldFields{ii});
    Subjects = rmfield(Subjects, oldFields{ii});
end
%%
start_prob_vec = [0.9, 0.95];

for subj_cnt = 1 : Total_Subject_Count
    for ss = 1 : 2
        start_prob = start_prob_vec(ss);
        ind_freq_start = find(Subjects(subj_cnt).relative_frequency.p_vec < start_prob, 1, 'last');
        if ~isempty(ind_freq_start)
            ind_freq_start = ind_freq_start +1;
        else
            ind_freq_start = length(all_stim_vec);
        end
        fieldname = sprintf('p%d_start', start_prob*100);
        Subjects(subj_cnt).(fieldname).amplitude = all_stim_vec(ind_freq_start);
        Subjects(subj_cnt).(fieldname).prob = Subjects(subj_cnt).relative_frequency.p_vec(ind_freq_start);
        Subjects(subj_cnt).(fieldname).response = virtstimulate(Subjects(subj_cnt).(fieldname).amplitude, Subjects(subj_cnt).subj_parameters, noisefree);
    end
end
%%
Subjects = rmfield(Subjects, 'new_start');
Subjects = orderfields(Subjects, {'subj_parameters', 'relative_frequency', 'MLE_0050uV', 'p90_start', 'p95_start', 'NL_50uV', 'NL_500uV'});

file_name = fullfile(data_folder, sprintf('Subjects_%d_reprocessed.mat', Total_Subject_Count));
save(file_name, 'Subjects', 'all_stim_vec', '-append');
