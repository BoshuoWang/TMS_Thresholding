function batch_reSubject(batch_id, num_subj_per_batch)
addpath('Statistical-MEP-Model');   % IO model path
addpath('Functions');               % Function and LUT path
if ~exist('DataStorage', 'dir')
    mkdir('DataStorage');           % Data storage path
end
if ~exist(fullfile('Figures', 'IO_curves'), 'dir')
    mkdir(fullfile('Figures', 'IO_curves'));
end
if ~exist('Total_Subject_Count', 'var')
    Total_Subject_Count = 25e3;    % 25k
end
data_file_name = sprintf('Subjects_%d.mat', Total_Subject_Count);
load( fullfile('DataStorage', data_file_name), 'Subjects', 'Total_Subject_Count', 'y_thresh');
%   Subject fields: subj_parameters, start_amplitude, thresh_x, MLE_0050V
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters:
stim_amp_eps = 0.001e-2;                            %
stim_amp_step = 0.2e-2;                             % 0.2% MSO
stim_amp_min = stim_amp_step;
stim_amp_max = 1.5;                                 % 150% MSO

multiplier = 5;                                                             % ##
all_trial_count = 10 * multiplier;
coarse_trial_count = 40 * multiplier;
fine_trial_count = 200 * multiplier;

all_stim_vec = (stim_amp_min : stim_amp_step : stim_amp_max);               % 0%-150% MSO, in steps of 0.2% MSO (750 pulse amplitude)
all_stim_mat = repmat(all_stim_vec, [all_trial_count, 1]);                  % 0%-150% MSO, 10 trials at each pulse amplitude (7500 trials)
coarse_stim_offset = (-15e-2:stim_amp_step:15e-2);                          % 15% MSO below and above previous threshold, additional 40 trials (151 pulse amplitudes, 6040 trials total)
fine_stim_offset = (-3e-2:stim_amp_step:3e-2);                              % 3% MSO below and above previous threshold, additional 200 trials ( 31 pulse amplitudes, 6400 trials total)
ind_fit_offset = (abs(fine_stim_offset) <= 2e-2+stim_amp_eps);

start_prob = 0.9;                                                           % ##

% fsolve_opts = optimset('Display', 'off', 'MaxIter', 10000, 'TolFun', 1E-7, 'TolX', 1E-9);
MLE_opts = optimset('Display', 'off', 'MaxFunEvals', 500000, 'FunValCheck', 'on', 'MaxIter', 10000, 'TolFun', 1e-8, 'TolX', 1e-12);
LB = [0.01, 0.005];   % lower bounds for threshold and slope in MLE
UB = [1.3, 0.5];      % upper bounds for threshold and slope in MLE

noisefree = true;
% if_save = true;

T_subj_end = NaN(1, num_subj_per_batch);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subj_id_vec = (1 : num_subj_per_batch) + (batch_id - 1) * num_subj_per_batch;

NewSubjects = Subjects(subj_id_vec);
% oldFields = {'thresh_x', 'MLE_0050uV'}; ...'start_amplitude',
% newFields = {'NL50uV_thresh_x', 'old_MLE_0050uV'}; ...'old_start_amplitude',
% for ii = 1 : length(oldFields)
%     [NewSubjects.(newFields{ii})] = NewSubjects.(oldFields{ii});
%     NewSubjects = rmfield(NewSubjects, oldFields{ii});
% end
MLE_field_name = sprintf('MLE_%04duV', fix(y_thresh/1e-6));

%%
for subj_cnt = num_subj_per_batch : -1 : 1        % Reverse loop, elimitnate need to pre-allocate space for many variables
    subj_id = subj_id_vec(subj_cnt);
    fprintf('Processing subject No. %d\n', subj_id_vec(subj_cnt));
    T_subj_start = tic;

    subj_parameters = NewSubjects(subj_cnt).subj_parameters;

    NL50uV_thresh_x = NewSubjects(subj_cnt).thresh_x;

    fine_stim_vec = round(NL50uV_thresh_x/stim_amp_step)*stim_amp_step + fine_stim_offset;
    fine_stim_mat = repmat(fine_stim_vec, fine_trial_count, 1);
    fine_response_bin_mat = logical(virtstimulate(fine_stim_mat, subj_parameters) >= y_thresh);
    fine_response_count_vec = sum(fine_response_bin_mat, 1);
    ind_fine_stim = (all_stim_vec > fine_stim_vec(1)-stim_amp_eps) & (all_stim_vec < fine_stim_vec(end)+stim_amp_eps);

    coarse_stim_vec = round(NL50uV_thresh_x/stim_amp_step)*stim_amp_step + coarse_stim_offset;
    coarse_stim_mat = repmat(coarse_stim_vec, coarse_trial_count, 1);
    coarse_response_bin_mat = logical(virtstimulate(coarse_stim_mat, subj_parameters) >= y_thresh);
    coarse_response_count_vec = sum(coarse_response_bin_mat, 1);
    ind_coarse_stim = (all_stim_vec > coarse_stim_vec(1)-stim_amp_eps) & (all_stim_vec < coarse_stim_vec(end)+stim_amp_eps);

    all_response_mat = virtstimulate(all_stim_mat, subj_parameters);
    all_response_bin_mat = logical(all_response_mat >= y_thresh);
    all_response_count_vec = sum(all_response_bin_mat, 1);
    all_response_count_vec(ind_coarse_stim) = all_response_count_vec(ind_coarse_stim) + coarse_response_count_vec;
    all_response_count_vec(ind_fine_stim) = all_response_count_vec(ind_fine_stim) + fine_response_count_vec;
    
    all_freq_vec = zeros(size(all_response_count_vec));
    all_freq_vec(~ind_coarse_stim) = all_response_count_vec(~ind_coarse_stim)/all_trial_count;
    all_freq_vec(ind_coarse_stim & ~ind_fine_stim) = all_response_count_vec(ind_coarse_stim & ~ind_fine_stim)/(all_trial_count + coarse_trial_count);
    all_freq_vec(ind_fine_stim) = all_response_count_vec(ind_fine_stim)/(all_trial_count + coarse_trial_count + fine_trial_count);

    fine_freq_vec = all_freq_vec(ind_fine_stim);
    ind_th = find(fine_freq_vec>=0.50,1,'first');
    NewSubjects(subj_cnt).relative_frequency.p50_raw = fine_stim_vec(ind_th);
    stats = regstats(fine_freq_vec(ind_fit_offset), fine_stim_vec(ind_fit_offset), 'linear');
    interc = stats.beta(1);
    slope  = stats.beta(2);
    fine_freq_fit_vec = stats.yhat;

    NewSubjects(subj_cnt).relative_frequency.p_vec = all_freq_vec;
    NewSubjects(subj_cnt).relative_frequency.p50_lin = (0.5-interc)/slope;
    NewSubjects(subj_cnt).relative_frequency.response = virtstimulate((0.5-interc)/slope, subj_parameters, noisefree);
    NewSubjects(subj_cnt).relative_frequency.slope  = slope;
    NewSubjects(subj_cnt).relative_frequency.s_MLE = 1/(slope*sqrt(2*pi));
    NewSubjects(subj_cnt).relative_frequency.rsquare = stats.rsquare;
    
    NewSubjects(subj_cnt).NL50uV.thresh_x = NL50uV_thresh_x;
    NewSubjects(subj_cnt).NL50uV.prob = interc + NL50uV_thresh_x * slope;

    NewSubjects(subj_cnt).old_start.amplitude = NewSubjects(subj_cnt).start_amplitude;
    ind_old_start = find(all_stim_vec >= NewSubjects(subj_cnt).old_start.amplitude, 1, 'first');
    if ~isempty(ind_old_start)
        NewSubjects(subj_cnt).old_start.prob = all_freq_vec(ind_old_start);
    else
        NewSubjects(subj_cnt).old_start.prob = all_freq_vec(end);
    end

    ind_freq_start = find(all_freq_vec<start_prob, 1, 'last');
    if ~isempty(ind_freq_start)
        ind_freq_start = ind_freq_start +1;
    else
        ind_freq_start = length(all_stim_vec);
    end
    NewSubjects(subj_cnt).new_start.amplitude = all_stim_vec(ind_freq_start);
    NewSubjects(subj_cnt).new_start.prob = all_freq_vec(ind_freq_start);
    NewSubjects(subj_cnt).new_start.response = virtstimulate(NewSubjects(subj_cnt).new_start.amplitude, subj_parameters, noisefree);

    fun_MLE =  @(theta) (-1)*loglikelyhood([all_stim_mat(:)', coarse_stim_mat(:)', fine_stim_mat(:)'], [all_response_bin_mat(:)', coarse_response_bin_mat(:)', fine_response_bin_mat(:)'], theta(1), theta(2)); % theta(1) = t; theta(2) = s;
    theta_init = [0.5, 0.2];
    [theta, MLE_fval, MLE_exitflag] = fmincon(fun_MLE, theta_init, [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
    %                                 fmincon(FUN,     X0,         A,   B,  Aeq,Beq,LB,           UB,       NONLCON, OPTIONS)
    if ~MLE_exitflag
        count_run = 1;
        fprintf('MLE not successful. Rerun No. %d...\n', count_run);
        theta_init =  [rand * 1, rand * 0.15];
        [theta, MLE_fval, MLE_exitflag] = fmincon(fun_MLE, theta_init, [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
        while ~MLE_exitflag
            count_run = count_run + 1;
            fprintf('MLE not successful again. Rerun No. %d...\n', count_run);
            theta_init = rand * [rand * 1, rand * 0.1];
            [theta, MLE_fval, MLE_exitflag] = fmincon(fun_MLE, theta_init, [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
        end
    end
    NewSubjects(subj_cnt).(MLE_field_name).t = theta(1);
    NewSubjects(subj_cnt).(MLE_field_name).s = theta(2);
    NewSubjects(subj_cnt).(MLE_field_name).response = virtstimulate(theta(1), subj_parameters, noisefree);
    NewSubjects(subj_cnt).(MLE_field_name).MLE_fval = MLE_fval;
    NewSubjects(subj_cnt).(MLE_field_name).exitflag = MLE_exitflag;
    %%
    if ~exist('h_f', 'var') || ~isvalid(h_f)
        h_f = figure('Name', 'Threshold comparison', 'Color', 'w', 'Position', [0, 0, 2000, 800]);
        h_a0 = axes(h_f, 'Position', [0.1, 0.875, 0.825, 0.05] );
        axis(h_a0,'off');

        h_a(1) = axes(h_f, 'Position', [0.1, 0.105, 0.4, 0.8] );
        set(h_a(1), 'YScale', 'log', 'XLim', [0, stim_amp_max+stim_amp_step]*100, 'YLim', 3*[1e-6,1e-2]); % X: 0% to 100% MSO; Y: 3 �V to 30 mV
        xlabel(h_a(1), 'Stimulation amplitude (\% MSO)', 'Interpreter', 'Latex', 'FontSize', 20);
        ylabel(h_a(1), 'EMG response amplitude (V)', 'Interpreter', 'Latex', 'FontSize', 20);
        %  set(h_a(1), 'YTick', 10.^(-5:-2), 'YTickLabel', {'10 $\mathrm{\mu V}$', '100 $\mathrm{\mu V}$', '1 $\mathrm{mV}$', '10 $\mathrm{mV}$'});
        h_a(2) = axes(h_f, 'Position', [0.575, 0.105, 0.4, 0.8] );
        set(h_a(2), 'YScale', 'lin', 'YLim', [0,100]); % X: 0% to 100% MSO; Y: 3 �V to 30 mV
        xlabel(h_a(2),'Stimulation amplitude (\% MSO)', 'Interpreter', 'Latex', 'FontSize', 20);
        ylabel(h_a(2),'Suprathreshold response probability', 'Interpreter', 'Latex', 'FontSize', 20);
        set(h_a, 'Box', 'On', 'TickLabelInterpreter', 'Latex', 'FontSize', 18,'LineWidth', 1, 'NextPlot', 'Add');
        for ii = 1 : length(h_a(2).YTick)
            h_a(2).YTickLabel{ii} = sprintf('$$%g\\%%$$',h_a(2).YTick(ii));
        end

        h_a_in(1) = axes(h_f, 'Position', [0.12, 0.58, 0.1, 0.2] );
        set(h_a_in(1), 'YScale', 'log', 'YLim', 5e-5*10.^([-2,2]/50)); % X: 0% to 100% MSO; Y: 3 �V to 30 mV
        h_a_in(2) = axes(h_f, 'Position', [0.82, 0.155, 0.125, 0.25] );
        set(h_a_in(2), 'YScale', 'lin', 'YLim', 50+[-10,10]/4); % X: 0% to 100% MSO; Y: 3 �V to 30 mV
        set(h_a_in, 'Box', 'On', 'TickLabelInterpreter', 'Latex', 'FontSize', 14,'LineWidth', 1, 'NextPlot', 'Add');
        for ii = 1 : length(h_a_in(2).YTick)
            h_a_in(2).YTickLabel{ii} = sprintf('$$%g\\%%$$',h_a_in(2).YTick(ii));
        end
    end
    title(h_a0, sprintf('Subject No. %05i', subj_id), 'FontSize', 24, 'Interpreter', 'Latex')

    cla(h_a(1));
    plot(h_a(1), xlim(h_a(1)), y_thresh*[1,1]*10,'k:');
    plot(h_a(1), xlim(h_a(1)), y_thresh*[1,1],'k--');
    h_M = plot(h_a(1), all_stim_mat(:)*100, all_response_mat(:), '.', 'MarkerSize', 4, 'Color', [1,1,1]*0.5);
    h_s = plot(h_a(1), -1,-1, '.', 'MarkerSize', 12, 'Color', [1,1,1]*0.5);

    all_response_vec = virtstimulate(all_stim_vec, subj_parameters, noisefree);

    h_l1 = plot(h_a(1), all_stim_vec*100, all_response_vec, '-k', 'LineWidth', 2);
    h_p1 = plot(h_a(1), NL50uV_thresh_x*100, y_thresh, 'ok', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 2);
    h_p2 = plot(h_a(1), NewSubjects(subj_cnt).(MLE_field_name).t*100, NewSubjects(subj_cnt).(MLE_field_name).response, '^r', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.5);
    h_p3 = plot(h_a(1), NewSubjects(subj_cnt).relative_frequency.p50_lin*100, NewSubjects(subj_cnt).relative_frequency.response, 'vb', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.5);
    h_p4 = plot(h_a(1), NewSubjects(subj_cnt).old_start.amplitude*100, y_thresh*10, '<k', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 2);
    h_p5 = plot(h_a(1), NewSubjects(subj_cnt).new_start.amplitude*100, NewSubjects(subj_cnt).new_start.response, '>b', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 2);

    h_l = legend(h_a(1), [h_s, h_l1], { 'Samples', 'Noise-free IO curve'});
    set(h_l, 'Location', 'NorthWest', 'FontSize', 16, 'Box', 'On','Interpreter', 'Latex', 'EdgeColor', [1,1,1]);
    %     h_l = legend(h_a(1), [h_s, h_l1, h_p1, h_p2, h_p3, h_p4, h_p5], ...
    %         { 'Samples', 'Noise-free IO curve', 'Threshold: noise-free IO curve', 'Threshold: MLE', 'Threshold: relative-frequency, linear fit',...
    %         'Old starting point at 500 $\mathrm{\mu}$V', ['New starting point at $\ge$',num2str(start_prob,'%g'),'\%']});
    %     set(h_l, 'Location', 'NorthWest', 'FontSize', 16, 'Box', 'On','Interpreter', 'Latex', 'EdgeColor', [1,1,1]);
    %
    cla(h_a_in(1));
    plot(h_a_in(1), xlim(h_a(1)), y_thresh*[1,1],'k--');
    h_l1 = plot(h_a_in(1), all_stim_vec*100, all_response_vec, '-k', 'LineWidth', 2);
    h_p1 = plot(h_a_in(1), NL50uV_thresh_x*100, y_thresh, 'ok', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 2);
    h_p2 = plot(h_a_in(1), NewSubjects(subj_cnt).(MLE_field_name).t*100, NewSubjects(subj_cnt).(MLE_field_name).response, '^r', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.5);
    h_p3 = plot(h_a_in(1), NewSubjects(subj_cnt).relative_frequency.p50_lin*100, NewSubjects(subj_cnt).relative_frequency.response, 'vb', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.5);
    xlim(h_a_in(1), (NL50uV_thresh_x+[-1,1]*(stim_amp_max/2/50))*100);


    cla(h_a(2));
    plot(h_a(2), xlim(h_a(1)), [50,50],'k--');
    %     xlim(h_a(2), [min(fine_stim_vec)-0.1, max(fine_stim_vec)+0.1]*100);

    h_l1=...
        plot(h_a(2), all_stim_vec(ind_fine_stim)*100, all_freq_vec(ind_fine_stim)*100, 'o', 'LineWidth', 0.5, 'MarkerSize', 5, 'Color', [1,1,1]*0.5, 'MarkerFaceColor', [1,1,1]*0.5);
    plot(h_a(2), all_stim_vec(ind_coarse_stim & ~ind_fine_stim)*100, all_freq_vec(ind_coarse_stim & ~ind_fine_stim)*100, 'o', 'LineWidth', 0.5, 'MarkerSize', 4, 'Color', [1,1,1]*0.5, 'MarkerFaceColor', [1,1,1]*0.5);
    plot(h_a(2), all_stim_vec(~ind_coarse_stim)*100, all_freq_vec(~ind_coarse_stim)*100, 'o', 'LineWidth', 0.5, 'MarkerSize', 2, 'Color', [1,1,1]*0.5, 'MarkerFaceColor', [1,1,1]*0.5);
    MLE_p_vec = normcdf(all_stim_vec, NewSubjects(subj_cnt).(MLE_field_name).t, NewSubjects(subj_cnt).(MLE_field_name).s);
    h_l3 = plot(h_a(2), all_stim_vec*100, MLE_p_vec*100, 'r-', 'LineWidth', 1.5, 'Color', [1,0,0]);
    h_l2 = plot(h_a(2), fine_stim_vec(ind_fit_offset)*100, fine_freq_fit_vec*100, 'b-', 'LineWidth', 1.5, 'Color', [0,0,1]);

    h_p3 = plot(h_a(2), NL50uV_thresh_x*100, NewSubjects(subj_cnt).NL50uV.prob*100, 'ok', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    h_p4 = plot(h_a(2), NewSubjects(subj_cnt).(MLE_field_name).t*100, 0.50*100, '^r', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    h_p2 = plot(h_a(2), NewSubjects(subj_cnt).relative_frequency.p50_lin*100, 0.5*100 , 'vb', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    h_p1 = plot(h_a(2), NewSubjects(subj_cnt).relative_frequency.p50_raw*100, fine_freq_vec(ind_th)*100 , '+b', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    h_p5 = plot(h_a(2), NewSubjects(subj_cnt).old_start.amplitude*100, NewSubjects(subj_cnt).old_start.prob*100, '<k', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.5);
    h_p6 = plot(h_a(2), NewSubjects(subj_cnt).new_start.amplitude*100, NewSubjects(subj_cnt).new_start.prob*100, '>b', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.5);

    h_l = legend(h_a(2), [h_l1, h_l2, h_l3, h_p3, h_p4, h_p2, h_p1, h_p6, h_p5],...
        {'Relative frequency', 'Relative frequency: linear fit', 'Relative frequency: MLE fit', ...
        'Threshold: noise-free IO curve', 'Threshold: MLE', 'Threshold: relative frequency, linear fit', 'Threshold: relative frequency, raw',...
        ['New starting point at $\ge$',num2str(start_prob,'%g'),'\%'], 'Old starting point at 500 $\mathrm{\mu}$V'});

    set(h_l, 'Location', 'NorthWest', 'FontSize', 16, 'Box', 'off','Interpreter', 'Latex', 'EdgeColor', [1,1,1]);
    xlim(h_a(2), [min(coarse_stim_vec), max(coarse_stim_vec)]*100);

    cla(h_a_in(2));
    plot(h_a_in(2), xlim(h_a(2)), [50,50],'k--');
    xlim(h_a_in(2), [NL50uV_thresh_x-0.03/4, NL50uV_thresh_x+0.03/4]*100);
    h_l3 = plot(h_a_in(2), all_stim_vec*100, MLE_p_vec*100, 'r-', 'LineWidth', 1.5, 'Color', [1,0,0]);
    h_l2 = plot(h_a_in(2), fine_stim_vec(ind_fit_offset)*100, fine_freq_fit_vec*100, 'b-', 'LineWidth', 1.5, 'Color', [0,0,1]);
    h_p3 = plot(h_a_in(2), NL50uV_thresh_x*100, NewSubjects(subj_cnt).NL50uV.prob*100, 'ok', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    h_p4 = plot(h_a_in(2), NewSubjects(subj_cnt).(MLE_field_name).t*100, 0.50*100, '^r', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    h_p2 = plot(h_a_in(2), NewSubjects(subj_cnt).relative_frequency.p50_lin*100, 0.5*100 , 'vb', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'w');
    h_p1 = plot(h_a_in(2), NewSubjects(subj_cnt).relative_frequency.p50_raw*100, fine_freq_vec(ind_th)*100 , '+b', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'w');

    figure_name = fullfile('Figures', 'IO_curves', sprintf('Subject_%05i', subj_id) );
    saveas(gcf, [figure_name, '.png']);

    T_subj_end(subj_cnt) = toc(T_subj_start);
    T_run_min = sum(T_subj_end, 'omitnan')/60;
    T_subj_end_ave = mean(T_subj_end, 'omitnan');
    T_rem_min_est = num_subj_per_batch * T_subj_end_ave/60 - T_run_min;
    fprintf('\tSubject %03d:%3d:%02d:%06.3f.\tRun time: %3d hours %04.1f minutes.\n', subj_id, ...
        floor(T_subj_end(subj_cnt)/3600), floor(mod(T_subj_end(subj_cnt)/60, 60)), mod(T_subj_end(subj_cnt), 60), ...
        floor(T_run_min/60), mod(T_run_min,60));
    fprintf('\tAverage: %3d:%02d:%06.3f.\n',...
        floor(T_subj_end_ave/3600),  floor(mod(T_subj_end_ave/60, 60)),  mod(T_subj_end_ave, 60));
    if subj_cnt > 1
        fprintf('\tRemaining subjects: %3d.\tRemaining time: %3d hours %04.1f minutes (estimated).\n', subj_cnt-1, ...
            floor(T_rem_min_est/60), mod(T_rem_min_est, 60) );
    end
end
NewSubjects = rmfield(NewSubjects, {'start_amplitude', 'thresh_x'});

data_file_name = sprintf('reSubjects_No%d.mat', batch_id);
folder_name = fullfile('DataStorage', sprintf('reSubjects_%dsubjperbatch', num_subj_per_batch) );
if ~exist(folder_name, 'dir')
    mkdir(folder_name);           % Data storage path
end
save( fullfile(folder_name, data_file_name), 'NewSubjects');

fprintf('Saved.\n');
end