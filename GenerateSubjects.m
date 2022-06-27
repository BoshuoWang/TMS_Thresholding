addpath('Statistical-MEP-Model');   % IO model path
addpath('Functions');               % Function and LUT path
if ~exist('DataStorage','dir')
    mkdir('DataStorage');           % Data storage path
end
if ~exist(fullfile('Figures','IO_curves'), 'dir')
    mkdir(fullfile('Figures','IO_curves'));
end

name_prefix = 'cluster_run';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parameters:
y_thresh = 50e-6;                       % 50 µV, in V
stim_vec = [eps, 1.9 * rand([1, 10000]) ];% 0%-190% MSO, same for every subject
stim_vec_lin= linspace(eps, 1.5, 10000);

fsolve_opts = optimset('Display', 'off', 'MaxIter', 10000, 'TolFun', 1E-6, 'TolX', 1E-8);
MLE_opts = optimset('Display', 'off', 'MaxFunEvals', 500000, 'FunValCheck', 'on', 'MaxIter', 10000, 'TolFun', 1e-6, 'TolX', 1e-8);
LB = [0.01,0.001];   % Lower bounds for threshold and slope in MLE
UB = [1.5,0.5];      % Upper bounds for threshold and slope in MLE

if ~exist('Total_Subject_Count', 'var')
    Total_Subject_Count = 25e3;    % 25k
end
response_lvl_list = [y_thresh, 500e-6];

Subjects = struct('subj_parameters',cell(1,Total_Subject_Count),'thresh_x',[], 'start_amplitude',[],...
                  'MLE_0050uV',struct('t',[],'s',[],'MLE_fval',[]','exitflag',[]));

noisefree = true;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Main Loop:
if isunix % Set up parpool
    numCPUs = str2double(getenv('SLURM_CPUS_PER_TASK'));
    fprintf('\nSetting up cluster pool.\nNumber of CPUs requested = %d.\n', numCPUs);
    
    pc_storage_dir = fullfile('PC_storage',getenv('SLURM_JOB_ID'));    % assign JobStorageLocation based on id of SLURM job that called this function
    mkdir(pc_storage_dir);
    
    pc = parcluster('local');
    pc.JobStorageLocation =  pc_storage_dir;
    poolobj = parpool(pc, numCPUs-1);
else
    pc = parcluster('local');
end

parfor subj_cnt = 1 : Total_Subject_Count
    
    figure('Name', 'IO Curve', 'Color', 'w', 'Units', 'Normalized', 'Position', [0.05, 0.05, 0.85, 0.85]);
    set(gca, 'Box', 'On', 'TickLabelInterpreter', 'Latex', 'FontSize', 16,'LineWidth', 1,'NextPlot', 'Add');
    set(gca, 'YScale', 'log', 'XLim', [0, 150], 'YLim', 3*[1e-6,1e-2]); % X: 0% to 100% MSO; Y: 3 µV to 30 mV
    xlabel('Stimulation amplitude (\% MSO)', 'Interpreter', 'Latex', 'FontSize', 18);
    ylabel('Response amplitude (V)', 'Interpreter', 'Latex', 'FontSize', 18);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Subject parameters and threshold
    
    %   Generate paramameters for subject
    subjected_generated = false;
    
    while (~subjected_generated)
        subj_parameters = virtualsubjectEIVGenerateSubject();
        Subjects(subj_cnt).subj_parameters = subj_parameters;    % Store them for later
        %   Get expected threshold amplitudes and threshold from estimator with
        %   10000 points or more for several levels: 50 µV, 500 µV
        resp_vec = virtstimulate(stim_vec, subj_parameters);
        for response_lvl_ind = 1 : length(response_lvl_list)
            response_lvl = response_lvl_list(response_lvl_ind);
            fun_opt = @(stim_amp) virtstimulate(stim_amp, subj_parameters, noisefree) - response_lvl;
            x_init = min( max(response_lvl/2e2, 0.7), 1);
            yx_0_init = - response_lvl  + 10e-6;
            [   tmp_x, tmp_y, exit_flag ] = fsolve_diffser(fun_opt, x_init, yx_0_init, fsolve_opts);
            switch response_lvl_ind
                case 1
                    Subjects(subj_cnt).thresh_x  = tmp_x;
                    response_bin = logical(resp_vec >= response_lvl);
                    fun_MLE =  @(theta) (-1)*loglikelyhood(stim_vec, response_bin, theta(1), theta(2)); % theta(1) = t; theta(2) = s;
                    theta_init = [0.5, 0.2];
                    [theta, MLE_fval, MLE_exitflag] = fmincon(fun_MLE, theta_init, [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
                    %                                 fmincon(FUN,     X0,         A,  B,  Aeq,Beq,LB, UB, NONLCON, OPTIONS)
                    if ~MLE_exitflag
                        count_run = 1;
                        fprintf('MLE not successful. Rerun No. %d...\n', count_run);
                        theta_init =  [rand * 1, rand * 0.15];
                        [theta, MLE_fval, MLE_exitflag] = fmincon(fun_MLE, theta_init, [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
                        while ~MLE_exitflag
                            count_run = count_run + 1
                            fprintf('MLE not successful again. Rerun No. %d...\n', count_run);
                            theta_init = rand * [rand * 1, rand * 0.1];
                            [theta, MLE_fval, MLE_exitflag] = fmincon(fun_MLE, theta_init, [], [], [], [], LB, UB, [], MLE_opts);  % find max via min of negative
                        end
                    end
                    field_name = sprintf('MLE_%04duV', fix(response_lvl/1e-6));
                    Subjects(subj_cnt).(field_name).t = theta(1);
                    Subjects(subj_cnt).(field_name).s = theta(2);
                    Subjects(subj_cnt).(field_name).MLE_fval = MLE_fval;
                    Subjects(subj_cnt).(field_name).exitflag = MLE_exitflag;
                case 2
                    Subjects(subj_cnt).start_amplitude  = tmp_x;   % Always above threshold so that IFCN 5/10 or 10/20 works
                    
            end
        end
        fprintf('\tSubject No. %d: Threshold %5.1f%% MSO. \n', subj_cnt, Subjects(subj_cnt).thresh_x * 100);
        
        if (Subjects(subj_cnt).thresh_x > 0.3 && Subjects(subj_cnt).thresh_x < 1)
            subjected_generated = true;
            
            cla
            plot(stim_vec*100, resp_vec, '.k', 'MarkerSize', 5)
            plot(stim_vec_lin*100, virtstimulate(stim_vec_lin, subj_parameters, noisefree), '-r', 'LineWidth', 2)
            plot(Subjects(subj_cnt).thresh_x*100, y_thresh , 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
            
            field_name = sprintf('MLE_%04duV',  fix(y_thresh/1e-6));
            plot(Subjects(subj_cnt).(field_name).t*100, y_thresh, 'db', 'MarkerFaceColor', 'b');
            
            h_l = legend({'Samples', 'Noise-free IO curve', 'Thresholds', 'Estimated Thresholds'});
            set(h_l, 'Location', 'NorthWest', 'FontSize', 16, 'Box', 'On','Interpreter', 'Latex', 'EdgeColor', [1,1,1]);
            title(sprintf('Subject No. %04i', subj_cnt), 'FontSize', 18, 'Interpreter', 'Latex')
            figure_name = fullfile('Figures','IO_curves', sprintf('Subject_%05i', subj_cnt) );
            saveas(gcf, [figure_name, '.png']);
        else
            fprintf('\tSubject No. %d: Threshold outside range from 35%% MSO to 100%% MSO. Regenerate subject parameters. \n', subj_cnt);
        end
    end
    
end

if exist('poolobj','var')
    delete(poolobj);
    rmdir(pc_storage_dir, 's'); % delete parfor temporary files
end

data_file_name = sprintf('Subjects_%d.mat', Total_Subject_Count);
save( fullfile('DataStorage', data_file_name), ...
        'Subjects', 'Total_Subject_Count', 'y_thresh' );
 
fprintf('Saved.\n');

rmpath('Statistical-MEP-Model');   % IO model path
rmpath('Functions');               % Function and LUT path