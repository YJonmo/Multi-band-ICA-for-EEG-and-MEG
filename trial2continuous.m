function [data_continuous] = trial2continuous(data_trial)
% This function converts the trialed data to continius data. 
data_continuous = data_trial;
if  length(data_trial.trial) == 1
    data_trial.trialinfo = 1
end
data_continuous.trial_original = data_continuous.trial;
data_continuous.time_original = data_continuous.time;
data_continuous = rmfield(data_continuous, 'trial');
data_continuous = rmfield(data_continuous, 'time');
%template_time = (0: (data_trial.time{1}(2) - data_trial.time{1}(1)) : 1/data_trial.fsample*size(data_trial.time{1},2)) ;
template_time = linspace(0, data_trial.time{1}(end), size(data_trial.time{1},2)) ; 
data_continuous.trial{1}  = zeros(size(data_trial.trial{1},1),  size(data_trial.trial{1},2) * size(data_trial.trial,2));
data_continuous.time{1}  = zeros(1,                             size(data_trial.time{1},2)  * size(data_trial.time ,2));
for Trial_Index = 1:size(data_trial.trial,2)
    data_continuous.trial{1}(:,(Trial_Index-1)*size(data_trial.trial{Trial_Index},2)+ 1: Trial_Index*size(data_trial.trial{Trial_Index},2)) = data_trial.trial{Trial_Index}(:,:); 
    
    if  isfield(data_trial, 'data_trial.trialinfo')
        data_continuous.time{1}(1,(Trial_Index-1)*size(data_trial.time{Trial_Index},2)+ 1: Trial_Index*size(data_trial.time{Trial_Index},2)) = template_time(1:end-1) + (data_trial.trialinfo(Trial_Index)-1)*template_time(end) ;
    else
        
        if  Trial_Index == 1
            data_continuous.time{1}(1, 1:size(data_trial.time{1},2)) = template_time(1:end) ;
        else
            data_continuous.time{1}(1,(Trial_Index-1)*size(data_trial.time{Trial_Index},2)+ 1: Trial_Index*size(data_trial.time{Trial_Index},2)) = data_continuous.time{1}(1, (Trial_Index-1)*size(data_trial.time{Trial_Index-1},2)) + template_time(1:end) ;
        end
    end
end
if  isfield(data_trial, 'sampleinfo')
    data_continuous = rmfield(data_continuous, 'sampleinfo');
    data_continuous.sampleinfo(1) = data_trial.sampleinfo(1,1);
    %data_continuous.sampleinfo(2) = data_trial.sampleinfo(end,2);
    data_continuous.sampleinfo(2) = size(data_continuous.trial{1},2) ;    
end

%data_continuous.sampleinfo(3) = data_trial.sampleinfo(end,3);