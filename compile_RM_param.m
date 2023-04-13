function compile_RM_param(version, is_analog, num_subj_per_batch, num_files)
data_folder = 'DataStorage';

fprintf('Retrieving RM param data for version %d from: %s\n', version, data_folder);

if nargin <= 2
    num_subj_per_batch = 100;
    num_files = 250;
end
step_number = 200;

folder_name = fullfile(data_folder,sprintf('RobMon_%dsubjperbatch', num_subj_per_batch) );
if is_analog
    var_name =  'RobMonA';
else
    var_name =  'RobMonD';
end
err_str ={'abs','rel'};
Tdata = tic;
data_name = sprintf('%s_%d', var_name, version);

if exist(fullfile(data_folder, data_name), 'dir') == 0
    mkdir(fullfile(data_folder, data_name));
end

%%
for Batch_id = num_files : -1 : 1
    fprintf('\tProcessing file No.%d.\n', num_files+1-Batch_id)
    ind = (1 : num_subj_per_batch) + (Batch_id -1) * num_subj_per_batch ;
    
    if ~(version == 3 && is_analog)
        data_file_name = sprintf('RobMon_v%d_No%d.mat', version, Batch_id);
    else
        data_file_name = sprintf('RobMon_v%d_No%d.mat', version+1, Batch_id);
    end
    Data = load( fullfile(folder_name, data_file_name), [var_name,'*']);
    %%
    if Batch_id == num_files
        field_name = sprintf('%s%d_s', var_name, version);
        param_num_s = 21;
        param_num_w = 1;
        ind_w = 1;
        switch version
            case {1, 5, 7, 9, 73, 11}                
            case {3}
                if is_analog
                    field_name = sprintf('%s%d_sw', var_name, version+1);
                    num_second_order_weights = 41;
                    ind_w = (num_second_order_weights+1)/2;
                end
            case {4}
                field_name = sprintf('%s%d_sw', var_name, version);
                param_num_w = 41;
                ind_w = 1 : param_num_w;
        end
        param_num = param_num_s * param_num_w;                
    
        all_sub_field_names = fieldnames(Data.(field_name));
        for err_ind = 1 : 2
            for jj = 1 : length(all_sub_field_names)
                if strcmp(all_sub_field_names{jj}, [err_str{err_ind},'_err'])
                    break
                end
            end
            err_sub_field_names{err_ind} = all_sub_field_names([1:jj-1,jj+1:end]);
        end
        save(fullfile(data_folder, sprintf('Compiled_%dsubjs_%s_params.mat', num_subj_per_batch*num_files,  data_name)),...
                'param_num', 'param_num_s','param_num_w');
    end
    %%
    for err_ind = 1:2
        tmp_data = rmfield(Data.(field_name),err_sub_field_names{err_ind});
        if ~iscolumn(Data.(field_name))
            tmp_data = tmp_data';
        end
        tmp_data = permute(tmp_data,[3,2,1]);
        tmp_data = cell2mat(struct2cell(tmp_data));
        tmp_data = permute(tmp_data,[2,3,4,1]);
        eval(sprintf('%s_err_all(1 : param_num_s, 1 : param_num_w, ind, 1 : step_number +1) = tmp_data(:, ind_w, :, :);', err_str{err_ind}));
    end
   
end
T = toc(Tdata);
fprintf('Time for loading all data: %f s.\n\n', T);
Tdata = tic;
%%
for err_ind = 1 : 2    
    for jj = param_num_s : -1 : 1
        for ll = param_num_w : -1 : 1
            
            eval(sprintf('err = sort(squeeze(%s_err_all(jj, ll,:,:)),1);', err_str{err_ind}));
            err_stat(jj, ll).quartiles = prctile(err,(25:25:75),1);
            err_stat(jj, ll).max_err = max(err,[],1);
            err_stat(jj, ll).min_err = min(err,[],1);
            err_stat(jj, ll).median_abs = median(abs(err),1);
            
            for ii =  step_number + 1 : -1 : 1
                IQL = err_stat(jj, ll).quartiles(3,ii) - err_stat(jj, ll).quartiles(1,ii);
                err_stat(jj, ll).upperwisker(ii) = err(find(err(:,ii) <= err_stat(jj, ll).quartiles(3,ii) + 1.5 * IQL, 1, 'last' ),ii);
                err_stat(jj, ll).lowerwisker(ii) = err(find(err(:,ii) >= err_stat(jj, ll).quartiles(1,ii) - 1.5 * IQL, 1, 'first'),ii);
                err_stat(jj, ll).upperoutlier(ii) = sum(err(:,ii) > err_stat(jj, ll).quartiles(3,ii) + 1.5 * IQL);
                err_stat(jj, ll).loweroutlier(ii) = sum(err(:,ii) < err_stat(jj, ll).quartiles(1,ii) - 1.5 * IQL);
            end
            eval(sprintf('%s_err = err;', err_str{err_ind}));
            
            file_name = fullfile(data_folder, data_name, sprintf('Compiled_%d_subjs_dataset_No%03d.mat', num_subj_per_batch*num_files, (jj-1)*param_num_w+ll));
            
            starting_param = jj;
            weight_param = ll;
            if exist(file_name, 'file')
                append_flag = '-append';
            else
                append_flag = '';
            end
            save(file_name, sprintf('%s_err', err_str{err_ind}), 'starting_param', 'weight_param', append_flag);
            fprintf('Saved %s error data for parameter No. %d.\n', err_str{err_ind}, (jj-1)*param_num_w+ll);
%    
        end
    end
    eval(sprintf('%s_err_stat = err_stat;', err_str{err_ind}));
end
T = toc(Tdata);
fprintf('Time for saving results of individual parameters: %f s.\n\n', T);
save(fullfile(data_folder, sprintf('Compiled_%dsubjs_%s_err_stat.mat', num_subj_per_batch*num_files, data_name)),...
    '*_err_stat', 'param_num', 'param_num_s', 'param_num_w');
fprintf('\tSaved error stat.\n\n');
   
end
