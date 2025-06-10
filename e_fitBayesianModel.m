clear; clc; 
addpath('C:/Users/nnu02/Documents/MATLAB/fieldtrip-20220827');

subjects = ["04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "24" "25"];
%subjects = ["04"];

%set up paths
folder = 'E:/Gian/GG_SensAtt_Prediction/02Data';
touch_events = [125 124 252];
touch_labels = ["noTouch" "Touch"];

move_events = [2, 4; 1, 3];
move_labels = ["Stay", "Move"];

prob_events = ["25/75", "50/50", "75/25"];
prob_labels = ["LowProb", "EqualProb", "HighProb"];

k_all = []; 
k_true = []; 
Correctness_all = []; 

for ID = subjects
    ID = char(ID);

    subj_folder = strcat(folder, '/ID', ID, '/01EEG/');
    bdf_file = strcat(folder, '/ID', ID, '/01EEG/*ID', ID, '*.bdf');
    bdf_file = dir(bdf_file); 
    bdf_file = strcat(bdf_file.folder, '\', bdf_file.name); 


    %% Import data from EEG and deleted trials
    cfg                         = [];
    cfg.dataset                 = bdf_file;
    cfg.trialfun                = 'ft_trialfun_general'; % this is the default
    cfg.trialdef.eventtype      = 'STATUS';
    cfg.trialdef.eventvalue     = touch_events; % the values of the stimulus trigger for the three conditions
    cfg.trialdef.prestim        = 3; % in seconds
    cfg.trialdef.poststim       = 3; % in seconds
    cfg = ft_definetrial(cfg);
    trl = cfg.trl;
    events_list_EEG = [cfg.event.value];
    
    load(strcat(subj_folder, 'preprocessed_final_SPM_vEOG.mat'));
    load(strcat(subj_folder, 'rejected_trials.mat'));
    load(strcat(subj_folder, 'rejected_trials_RT'));


    %% Import Log
    opts = delimitedTextImportOptions("NumVariables", 7);
    opts.DataLines = [4, Inf];
    opts.Delimiter = ["\t", " \t "];
    opts.VariableNames = ["Block", "Trial", "BlockType", "TrialType", "Stimulation", "Event_Name", "Time"];
    opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts = setvaropts(opts, ["Block", "Trial", "BlockType", "TrialType", "Stimulation", "Event_Name", "Time"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Block", "Trial", "BlockType", "TrialType", "Stimulation", "Event_Name", "Time"], "EmptyFieldRule", "auto");
    Log = readmatrix(strcat("E:\Gian\GG_SensAtt_Prediction\02Data\ID", ID, "\00Behavioural\Experiment_Log_ID", ID, ".txt"), opts);
    for i = 1:length(Log)
        if contains(Log(i,1), "Start")
            Log(i,6) = Log(i,1);
            Log(i,7) = Log(i,2);
            Log(i,1) = Log(i+1,1);
            Log(i,2) = Log(i+1,2);
            Log(i,3) = Log(i+1,3);
            Log(i,4) = Log(i+1,4);
            Log(i,5) = Log(i+1,5);
        end
    end

    %first "correct" the event file with the Log file. This
    %makes sure that all events are in the correct order :)
    
    events_list_EEG(events_list_EEG == 126) = [];
    events_list_EEG(events_list_EEG~=10 & events_list_EEG~=20 & events_list_EEG~=30 & events_list_EEG~=124 & events_list_EEG~=125) = []; 
    events_list_EEG = events_list_EEG'; 

    events_log_list = Log(:,6);
    index = startsWith(events_log_list,"Start_trial");
    events_log_list(index) = "10";
    events_log_list = strrep(events_log_list,'preGO','20');
    events_log_list = strrep(events_log_list,'GO','30');
    events_log_list = strrep(events_log_list,'Index_visual','125');
    events_log_list = str2double(events_log_list);
    events_log_list(Log(:,6) == "Index_visual" & Log(:,5) == "1",:) = 124; 

    %check for repeated ones in both files
    repeated_EEG = []; 
    for i = 2:length(events_list_EEG)
        if events_list_EEG(i) == events_list_EEG(i-1)
            repeated_EEG = [repeated_EEG i]; 
        end
    end
    events_list_EEG(repeated_EEG) = []; 

    repeated_LOG = []; 
    for i = 2:length(events_log_list)
        if events_log_list(i) == events_log_list(i-1)
            repeated_LOG = [repeated_LOG i]; 
        end
    end
    events_log_list(repeated_LOG) = []; 
    Log1 = Log; 
    Log1(repeated_LOG,:) = []; 

    i = 1;
    errors = [];
    events_list_new = events_list_EEG;
    while (i < length(events_log_list)) 
        if events_list_new(i) == events_log_list(i)
            i = i + 1;
        else
            idx = i;
            events_list_new = [events_list_new(1:length(events_list_new) < idx); events_log_list(i); events_list_new(1:length(events_list_new) >= idx)];
            errors = [errors i];
            i = i+2;
        end
    end
    if sum(events_list_new ~= events_log_list)>0 
        error('something went terribly wrong');
    end

    Log1(errors,:) = []; 
    list_conditions_nonCorr = Log1(Log1(:,6) == "Index_visual",[4 5 2 1 3]);

    list_conditions = list_conditions_nonCorr; 
    list_conditions(rejected_trials_RT,:) = []; 

    %check for correctness of data
    zeros_ones_EEG = abs(data.trialinfo - 125);
    zeros_ones_Log = list_conditions_nonCorr(:,2); 
    zeros_ones_Log(unique([rejected_trials rejected_trials_RT])) = []; 
    zerosss = [zeros_ones_EEG, zeros_ones_Log];
    if sum(zeros_ones_EEG ~= str2double(zeros_ones_Log)) > 0
        error('PORCAMADONNA IL LOG E` SBAGLIATO')
    else
        display('All good!')
    end

    %% Import participant's responses
    opts = delimitedTextImportOptions("NumVariables", 2);
    opts.DataLines = [3, Inf];
    opts.Delimiter = "\t";
    opts.VariableNames = ["Var1", "Time"];
    opts.SelectedVariableNames = "Time";
    opts.VariableTypes = ["string", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");
    Experimentresponse = readtable(strcat("E:\Gian\GG_SensAtt_Prediction\02Data\ID", ID, "\00Behavioural\Experiment_response_ID", ID, ".txt"), opts);
    Experimentresponse = table2array(Experimentresponse);
    clear opts
    
    %import TrialTable from experiment
    opts = delimitedTextImportOptions("NumVariables", 13);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["Block", "Block_num", "TrialInBlock", "Trial", "Block_type", "Trial_type", "Stimulation", "ITIs", "Speed", "Completed", "Attempts", "Skipped", "TrialTime"];
    opts.VariableTypes = ["double", "double", "double", "double", "categorical", "double", "double", "double", "double", "categorical", "double", "categorical", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts = setvaropts(opts, ["Block_type", "Completed", "Skipped"], "EmptyFieldRule", "auto");
    TrialTable = readtable(strcat("E:\Gian\GG_SensAtt_Prediction\02Data\ID", ID, "\00Behavioural\ID", ID, "TrialTable.csv"), opts);
    clear opts

    ExperimentCondition = string(TrialTable{TrialTable.TrialInBlock == 0, 5});
    ExperimentCondition(ExperimentCondition == "25/75") = "25"; 
    ExperimentCondition(ExperimentCondition == "50/50") = "50"; 
    ExperimentCondition(ExperimentCondition == "75/25") = "75"; 
    ExperimentCondition = str2double(ExperimentCondition); 

    Correctness = ExperimentCondition == Experimentresponse; 
    Correctness_all = [Correctness_all sum(Correctness)/24];

    %% Adjust values for bayesian model fitting
    %first observations
    a = [];
    for i = 1:24
        a = [a sum(double(list_conditions(list_conditions(:,4) == string(i-1),2)))];
    end

    x = []; 
    x = Experimentresponse'/100;

    k0 = 1;

    k = fminsearch(@(k)negativeLog_BetaBernoulli(k,x,a),k0);

    x_real = ExperimentCondition'/100; 
    k_real = fminsearch(@(k)negativeLog_BetaBernoulli(k,x_real,a),k0);

    k_all = [k_all k];
    k_true = [k_true k_real]; 

end

%% PLOT STUFF
%first plot an example of beta updating
observations = [0 1 0 0 1 1 1 0 0 1 0 1 0 1 1 0 0 0 0 1 1 0 1 0 1];
x = linspace(0,1,1000);
figure;
for obs = 1:length(observations)
    color = [1-obs/length(observations), 1-obs/length(observations), 1];
    a = 1 + sum(observations(1:obs)==1);
    b = 1 + sum(observations(1:obs) == 0);
    plot(x, betapdf(x,a,b), 'Color', color)
    hold on
end
title('Updating of the belief over 25 observations, k = 1')
figure;
plot(observations,[1:25], 'b*')
xlim([-1 2]);
set(gca, 'YDir','reverse')


fplot(@(k)negativeLog_BetaBernoulli(k,x,a))

x = linspace(0,1,1000);
figure;
hold on
for k_example = 0.1:0.1:1
    color = [1-k_example 1-k_example 1];
    plot(x, betapdf(x, k_example*13, k_example*12), 'Color', color)
end
figure; 
k = 1;
a = 12;
b = 13;
plot(x, betapdf(x, k*a, k*b), 'Color', [0 0 1])
title(strcat('Beta(', 'k*(1+', int2str(a), ',),', 'k*(1+', int2str(b), ')); k = ', num2str(k)))
fplot(@(k)negativeLog_BetaBernoulli(k,x(1),a(1)))

%plot true values
figure; 
bar(subjects, k_true)
xlabel('Participants')
ylabel('k')
title('Barplot of the k values if participants answered correctly at all trials')

%plot correlation between correctness and K
mdl = fitlm(Correctness_all, k_all);
figure;
%plot(Correctness_all, k_all, 'k*')
plot(mdl)
ylim([0 2])
title(sprintf('Correlation between k and accuracy levels, p = %f', mdl.Coefficients{2,4}))
xlabel('Accuracy')
ylabel('k')
