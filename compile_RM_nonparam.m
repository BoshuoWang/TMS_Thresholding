function compile_RM_nonparam(version, num_subj_per_batch, num_files)
data_folder = 'DataStorage';

fprintf('Retrieving RM data from: %s\n',data_folder);

if nargin == 1
    num_subj_per_batch = 100;
    num_files = 250;
end
step_number = 200;

folder_name = fullfile('DataStorage',sprintf('RobMon_%dsubjperbatch', num_subj_per_batch) );
var_name =  'RobMon*';


for Batch_id = num_files : -1 : 1
    ind = (1 : num_subj_per_batch) + (Batch_id -1) * num_subj_per_batch ;

    data_file_name = sprintf('RobMon_v%d_No%d.mat', version, Batch_id);
    Data = load( fullfile(folder_name, data_file_name), var_name);


    switch version
        case {1, 5, 7, 9, 73, 11}
            num_start_conditions = 21;
            ind_start_conditions = (num_start_conditions+1)/2;
            eval(sprintf('RobMonA%d_abs_err(ind, 1 : step_number +1) = cell2mat(cellfun(@(x) x(:,ind_start_conditions), {Data.RobMonA%d_s.abs_err},''UniformOutput'',false))'';', version, version));
            eval(sprintf('RobMonA%d_rel_err(ind, 1 : step_number +1) = cell2mat(cellfun(@(x) x(:,ind_start_conditions), {Data.RobMonA%d_s.rel_err},''UniformOutput'',false))'';', version, version));
            eval(sprintf('RobMonA%d_run_time(ind, 1 : step_number) = cell2mat({Data.RobMonA%d_s.run_time})'';', version, version));
            if (version < 10)
                eval(sprintf('RobMonD%d_abs_err(ind, 1 : step_number +1) = cell2mat(cellfun(@(x) x(:,ind_start_conditions), {Data.RobMonD%d_s.abs_err},''UniformOutput'',false))'';', version, version));
                eval(sprintf('RobMonD%d_rel_err(ind, 1 : step_number +1) = cell2mat(cellfun(@(x) x(:,ind_start_conditions), {Data.RobMonD%d_s.rel_err},''UniformOutput'',false))'';', version, version));
                eval(sprintf('RobMonD%d_run_time(ind, 1 : step_number) = cell2mat({Data.RobMonD%d_s.run_time})'';', version, version));
            end
        case {3}
            eval(sprintf('RobMonA%d_MLE_abs_err(ind, 1 : step_number)  = cell2mat({Data.RobMonA%d.MLE_abs_err})'';', version, version));
            eval(sprintf('RobMonA%d_MLE_rel_err(ind, 1 : step_number)  = cell2mat({Data.RobMonA%d.MLE_rel_err})'';', version, version));
            eval(sprintf('RobMonA%d_lin_abs_err(ind, 1 : step_number)  = cell2mat({Data.RobMonA%d.lin_abs_err})'';', version, version));
            eval(sprintf('RobMonA%d_lin_rel_err(ind, 1 : step_number)  = cell2mat({Data.RobMonA%d.lin_rel_err})'';', version, version));

            num_start_conditions = 21;
            ind_start_conditions = (num_start_conditions+1)/2;
            eval(sprintf('RobMonD%d_abs_err(ind, 1 : step_number +1) = cell2mat(cellfun(@(x) x(:,ind_start_conditions), {Data.RobMonD%d_s.abs_err},''UniformOutput'',false))'';', version, version));
            eval(sprintf('RobMonD%d_rel_err(ind, 1 : step_number +1) = cell2mat(cellfun(@(x) x(:,ind_start_conditions), {Data.RobMonD%d_s.rel_err},''UniformOutput'',false))'';', version, version));
            eval(sprintf('RobMonD%d_run_time(ind, 1 : step_number) = cell2mat({Data.RobMonD%d_s.run_time})'';', version, version));

            data_file_name = sprintf('RobMon_v%d_No%d.mat', version+1, Batch_id);
            Data = load( fullfile(folder_name, data_file_name), var_name);
            num_second_order_weights = 41;
            ind_second_order_weights = (num_second_order_weights+1)/2;
            RobMonA3_abs_err(ind, 1 : step_number +1) = cell2mat(cellfun(@(x) x(:,ind_start_conditions, ind_second_order_weights),cellfun(@squeeze,{Data.RobMonA4_sw.abs_err},'UniformOutput',false),'UniformOutput',false))';
            RobMonA3_rel_err(ind, 1 : step_number +1) = cell2mat(cellfun(@(x) x(:,ind_start_conditions, ind_second_order_weights),cellfun(@squeeze,{Data.RobMonA4_sw.rel_err},'UniformOutput',false),'UniformOutput',false))';
            RobMonA3_run_time(ind, 1 : step_number) = cell2mat({Data.RobMonA4_sw.run_time})';

        case {2,6,8}
            eval(sprintf('RobMonA%d_abs_err(ind, 1 : step_number +1)  = cell2mat({Data.RobMonA%d.abs_err})'';', version, version));
            eval(sprintf('RobMonA%d_rel_err(ind, 1 : step_number +1)  = cell2mat({Data.RobMonA%d.rel_err})'';', version, version));
            eval(sprintf('RobMonA%d_run_time(ind, 1 : step_number) = cell2mat({Data.RobMonA%d.run_time})'';', version, version));
        case {4}
            num_start_conditions = 21;
            ind_start_conditions = (num_start_conditions+1)/2;
            num_second_order_weights = 41;
            default_second_order_weight = -0.10;
            second_order_weights = linspace(-1, 1, num_second_order_weights);
            ind_second_order_weights = find(abs(second_order_weights - default_second_order_weight)<eps);
            RobMonA4_abs_err(ind, 1 : step_number +1) = cell2mat(cellfun(@(x) x(:,ind_start_conditions, ind_second_order_weights),cellfun(@squeeze,{Data.RobMonA4_sw.abs_err},'UniformOutput',false),'UniformOutput',false))';
            RobMonA4_rel_err(ind, 1 : step_number +1) = cell2mat(cellfun(@(x) x(:,ind_start_conditions, ind_second_order_weights),cellfun(@squeeze,{Data.RobMonA4_sw.rel_err},'UniformOutput',false),'UniformOutput',false))';
            RobMonA4_run_time(ind, 1 : step_number) = cell2mat({Data.RobMonA4_sw.run_time})';
        case {13}
            RobMonA13_abs_err(ind, 1 : step_number +1) = cell2mat({Data.RobMonA13.abs_err})';
            RobMonA13_rel_err(ind, 1 : step_number +1) = cell2mat({Data.RobMonA13.rel_err})';
            RobMonA13_lin_abs_err(ind, 1 : step_number) = cell2mat({Data.RobMonA13.lin_abs_err})';
            RobMonA13_lin_rel_err(ind, 1 : step_number) = cell2mat({Data.RobMonA13.lin_rel_err})';
            RobMonA13_MLE_abs_err(ind, 1 : step_number) = cell2mat({Data.RobMonA13.MLE_abs_err})';
            RobMonA13_MLE_rel_err(ind, 1 : step_number) = cell2mat({Data.RobMonA13.MLE_rel_err})';
            RobMonA13_run_time(ind, 1 : step_number) = cell2mat({Data.RobMonA13.run_time})';
    end
end

file_name = fullfile(data_folder,sprintf('Compiled_%dsubjs_RM_%d.mat', num_subj_per_batch*num_files, version));
save(file_name,  var_name);
fprintf('Saved data to: %s\n',file_name);

end
