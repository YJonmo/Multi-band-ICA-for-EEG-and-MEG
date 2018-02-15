function MBICA = MultiBand_ICA_Simplified(cfg, dataToRun)

%% Minor manipulations and check ups
cfg.elec.label = lower(cfg.elec.label);
dataToRun.label = lower(dataToRun.label) ;   

if isfield(cfg, 'bandpass')
    if isfield(cfg, 'bandRanks')
        if size(cfg.bandpass,1) ~= size(cfg.bandRanks,1) 
            error('Error: The lenghts of bandpass and bandRanks should be the same')
        end
    end
    if isfield(cfg, 'rankThresh')
         if size(cfg.bandpass,1) ~= size(cfg.rankThresh,1) 
            error('Error: The lenghts of bandpass and rankThresh should be the same')
         end
    end
    if isfield(cfg, 'reducerankby')
         if size(cfg.bandpass,1) ~= size(cfg.rankThresh,1) 
            error('Error: The lenghts of bandpass and reducerankby should be the same')
         end
    end
end


%% Data cannot be too long for the sensor-space ICA (or it can contain all the trials but down sampling is needed)
% This part may reduce the EEG length as provided by the user.

if  isfield(cfg, 'trials')                                                % Reducing the duration of the data using predefined number of the trials.                   
    trials = cfg.trials;
    data_short = dataToRun ;
    data_short = rmfield(data_short, 'trial') ;
    data_short = rmfield(data_short, 'time') ;
    if isfield(dataToRun, 'trialinfo')    
        data_short = rmfield(data_short, 'trialinfo') ;
    end
    if isfield(dataToRun, 'sampleinfo')
        data_short = rmfield(data_short, 'sampleinfo') ;
    end
    for Trial_Index = trials(1):trials(end)
        data_short.trial{Trial_Index} = dataToRun.trial{Trial_Index};
        data_short.time{Trial_Index} = dataToRun.time{Trial_Index}; 

        if isfield(dataToRun, 'trialinfo')
            data_short.trialinfo(Trial_Index) = dataToRun.trialinfo(Trial_Index); 
        end
        if isfield(dataToRun, 'sampleinfo')
            data_short.sampleinfo(Trial_Index,:) = dataToRun.sampleinfo(Trial_Index,:);
        end
    end
    
    Sensor_space = trial2continuous(data_short);
    Sensor_space.continuous = Sensor_space.trial ;
            
elseif isfield(cfg, 'timewindow')                                        % Reducing the duration of the data using predefined number of seconds (start and end).
    Sensor_space = trial2continuous(dataToRun);
    Sensor_space.continuous = Sensor_space.trial ;
    Sensor_space.continuous{1} = Sensor_space.continuous{1}(:, 1 + cfg.timewindow(1)*dataToRun.fsample :  cfg.timewindow(2)*dataToRun.fsample) ;
    Sensor_space.time{1} = Sensor_space.time{1}(:, 1 + cfg.timewindow(1)*dataToRun.fsample :  cfg.timewindow(2)*dataToRun.fsample) ;
    Sensor_space.trial = Sensor_space.continuous;
    data_short = dataToRun ;
    data_short.trial = Sensor_space.continuous ; 
    data_short.trialinfo = 1 ;  % This means data is signle trial.
    data_short.time = Sensor_space.time ; 
    data_short.sampleinfo = [cfg.timewindow(1)*dataToRun.fsample+1 cfg.timewindow(2)*dataToRun.fsample] ; 
else
    data_short = dataToRun ;
    if  length(dataToRun.trial) == 1 ;
        data_short.trialinfo = 1 ;  % This means data is signle trial.
    elseif ~isfield(dataToRun, 'trial') 
        error('the input data is missing the .trial')
    end
    
    Sensor_space = trial2continuous(data_short);
    Sensor_space.continuous = Sensor_space.trial ;
end

No_Chan = size(dataToRun.trial{1},1) ;

%% Band-pass filtering the sensor space data
if  isfield(cfg, 'bandpass') 
    ft_progress('init', 'text', 'Band-pass filtering the sensor-space data ...') ; 
    if  isempty(cfg.bandpass)                                                                       
        Freqs =[0.5 4; 4 8; 8 12; 12 25; 25 47; 53 100]  ;
        %Names= { 'DeltaOut' 'ThetaOut' 'AlphaOut' 'BetaOut' 'BetaHighOut'}
    else
        Freqs = cfg.bandpass ; 
    end
    cfg_filt = [] ;
    cfg_filt.bpfilter = 'yes';
    Sensor_space_BPFiltered = Sensor_space ;
    Sensor_space_temp = Sensor_space;
    %Sensor_space_temp.time = Sensor_space.time_original ; 
    rmfield(Sensor_space_BPFiltered, 'continuous') ; 
    rmfield(Sensor_space_BPFiltered, 'trial') ; 
    rmfield(Sensor_space_temp, 'continuous') ; 
    rmfield(Sensor_space_temp, 'trial') ; 
    Sensor_space_temp.trial{1} = Sensor_space.continuous{1} ; 
    Sensor_space_BPFiltered.trial{1} = zeros(size(Sensor_space.continuous{1},1)*size(Freqs,1), size(Sensor_space.continuous{1},2)) ; 
    for Freq_Index = 1:size(Freqs,1)
        Freq = Freqs(Freq_Index,:) ; 
        ft_progress(Freq_Index/size(Freqs,1), 'Running the Band pass filter for frequency bands %d to %d', Freq(1), Freq(2));
        %Filter the data into the frrquency band of interest
        cfg_filt.bpfreq = Freq ;
        try
            [Sensor_space_temp2] = ft_preprocessing(cfg_filt, Sensor_space_temp) ;
        catch 
            if Freq_Index == 1;
                ft_progress('init', 'text', 'It is using the low pass for the first band rather than the defined band-pass.')
                cfg_filt = [] ;
                cfg_filt.lpfilter = 'yes';
                cfg_filt.lpfreq = Freqs(1,2)
                [Sensor_space_temp2] = ft_preprocessing(cfg_filt, Sensor_space_temp) ;    
                cfg_filt = [] ;
                cfg_filt.bpfilter = 'yes';
            end
        end
        Sensor_space_BPFiltered.Freqs{Freq_Index} = Sensor_space_temp2.trial{1} ; 
    end   
    Sensor_space = Sensor_space_BPFiltered ;    
    ft_progress('close') ; 
else 
    Sensor_space.Freqs{1} = Sensor_space.continuous{1} ; 
    Freq_Index = 1 ; 
end
rmfield(Sensor_space, 'continuous') ; 
clear Sensor_space_temp ; 
%% Runing the PCA

data_TemporalSubSpace = Sensor_space;
data_TemporalSubSpace = rmfield(data_TemporalSubSpace, 'trial') ;
data_TemporalSubSpace = rmfield(data_TemporalSubSpace, 'label') ;

if  isfield(cfg, 'bandpass') || isfield(cfg, 'bandRanks')           
    if isfield(cfg, 'bandRanks')                                        
        rank_source = cfg.bandRanks ;                               % User has defined the number of the PCs for each band
    else
        rank_source = zeros(1,size(Freqs,1)) - 1 ;                  % User has not defined the rank of each band and expects the automatic rank assignment
    end
else   
    ft_progress('init', 'text', 'Calculating the rank of the sensor-space data');
    rank_source = rank(Sensor_space.Freqs{1}) - 1                        % User has not defined the rank of each band but expect the rank for each bad to be reduced after automatic rank test
end

figure 
hold on

if isfield(cfg, 'rankThreshTotl')
    ft_progress('init', 'text','Calculating the rank of the sensor-space data, first frequency band')
    if isfield(cfg, 'pca')
        if sum(cfg.pca) == 333
            [U Sig V] = svd(Sensor_space.continuous{1}) ; 
            Eigens = diag(Sig).^2 ;
        elseif sum(cfg.pca) == 530
            cfg_temp.method = 'svd' ; 
            data_temp = [] ; 
            data_temp.label = Sensor_space.label ;
            data_temp.time{1} = Sensor_space.time{1} ;
            data_temp.trial{1} = Sensor_space.continuous{1} ;
            data_temp.fsample = Sensor_space.fsample ;
            [SVDed_data] = ft_componentanalysis(cfg_temp, data_temp) ; 
            Eigens = rms(SVDed_data.trial{1}') ;
            V = SVDed_data.trial{1}' ;
            U = SVDed_data.topo ; 
            Eigens = rms(SVDed_data.trial{1}') ;
            Eigen_Max = Eigens(1) ; 
        else
            [V U Eigens] = pca(Sensor_space.continuous{1}) ; 
        end            
    else
        [V U Eigens] = pca(Sensor_space.continuous{1}) ; 
    end
    Eigen_Max = Eigens(1) ; 
end
for SVD_Index = 1:Freq_Index      
%    ft_progress(SVD_Index/Freq_Index, 'Running the PCA of sensor-space data for frequency band %d from %d', SVD_Index, Freq_Index);
    if isfield(cfg, 'pca')
        if sum(cfg.pca) == 333
            [U Sig V] = svd(Sensor_space.Freqs{SVD_Index}) ; 
            Eigens = diag(Sig).^2 ;
        elseif sum(cfg.pca) == 646
            cfg_temp.method = 'svd' ; 
            data_temp = [] ; 
            data_temp.label = Sensor_space.label ;
            data_temp.time{1} = Sensor_space.time{1} ;
            data_temp.trial{1} = Sensor_space.Freqs{SVD_Index} ;
            data_temp.fsample = Sensor_space.fsample ;
            [SVDed_data] = ft_componentanalysis(cfg_temp, data_temp) ; 
            Eigens = rms(SVDed_data.trial{1}') ;
            V = SVDed_data.trial{1}' ;
            U = SVDed_data.topo ; 
            Eigens = rms(SVDed_data.trial{1}') ;
            Eigen_Max = Eigens(1) ; 
        else
            [V U Eigens] = pca(Sensor_space.Freqs{SVD_Index}) ; 
        end            
    else
        [V U Eigens] = pca(Sensor_space.Freqs{SVD_Index}) ; 
    end
    
    if  rank_source(SVD_Index) == -1
        %ft_progress('init', 'text', 'Calculating the rank of the sensor-space data');
        rank_source(SVD_Index) = rank(Sensor_space.Freqs{SVD_Index})   ;
        if isfield(cfg, 'rankThresh')
            for Eigen_index = 1:rank_source(SVD_Index)                              % Using the threshold in precentage of the first eigen value of each band to decide the rank of the source
                if  (Eigens(1,1)*cfg.rankThresh(SVD_Index))/100 > Eigens(Eigen_index,Eigen_index)
                    rank_source(SVD_Index) =  Eigen_index - 1 
                    break
                end
            end
        elseif isfield(cfg, 'rankThreshTotl')
            for Eigen_index = 1:rank_source(SVD_Index)                              % Using the threshold in precentage of the first eigen value of the wide-band to decide the rank of the source
                if  (Eigen_Max*cfg.rankThreshTotl)/100 > Sig(Eigen_index,Eigen_index)
                    rank_source(SVD_Index) =  Eigen_index - 1 
                    break
                end
            end      
        end
    end
    plot(Eigens)
    if  SVD_Index == 1
        Eigen_D = Eigens(1:rank_source(1)) ; 

        SpatialPCs(:, 1 : rank_source(1))  = U(:,1:rank_source(1)) ;       % Spatial subspace
        clear U ;
        TemporalSubSpace(:,:)  = V(:,1:rank_source(1)) ;                   % Temporal subspace  
        clear V ;
        clear Eigens
        for I = 1:rank_source(1)
            TemporalPCs(I,:) = sqrt(Eigen_D(I))*TemporalSubSpace(:,I)' ;   % Principal components (time-courses)
        end     
        data_TemporalSubSpace.trial{1}(1 : rank_source(1), :) = TemporalPCs(1 : rank_source(1),:) ;       
    else
        Eigen_D = Eigens(1:rank_source(SVD_Index)) ;       
        
        SpatialPCs(:, 1 + sum(rank_source(1:SVD_Index-1)): sum(rank_source(1:SVD_Index)))  = U(:,1:rank_source(SVD_Index)) ;      % Spatial subspace
        clear U ;
        TemporalSubSpace  = V(:,1:rank_source(SVD_Index)) ;                % Temporal subspace  
        clear V ;
        clear Eigens
        I2 = 1 ;
        for I = 1 + sum(rank_source(1:SVD_Index-1)): sum(rank_source(1:SVD_Index))
            TemporalPCs(I,:) = sqrt(Eigen_D(I2))*TemporalSubSpace(:,I2)' ; % principal components (time-courses)
            %TemporalPCs(I,:) = TemporalSubSpace(:,I2)' ; 
            I2 = I2 + 1 ;
        end        
        data_TemporalSubSpace.trial{1}(1 + sum(rank_source(1:SVD_Index-1)): sum(rank_source(1:SVD_Index)), :) = TemporalPCs(1 + sum(rank_source(1:SVD_Index-1)): sum(rank_source(1:SVD_Index)),:) ;
    end
end
%ft_progress('close') ;
for I = 1:Freq_Index
    plotlabel{I} = strcat('Band', num2str(I)) ; 
end
legend(plotlabel)
title('Eigen values of different bands')

%for I = 1:rank_source*Freq_Index
for I = 1:sum(rank_source)
    data_TemporalSubSpace.label{I} = strcat('PC',num2str(I));
end

%% Running the ICA

       
cfg2 = [];
TICA = data_TemporalSubSpace ; 
if  isfield(cfg, 'method')
    cfg2.method = cfg.method ;
    if  sum(cfg.method) == sum('sobi')
        if  isfield(cfg, 'sobi_delay')
            sobi_delay = cfg.sobi_delay ; 
        else 
            sobi_delay = 1 ; 
        end
            [Mixing trial] = sobi(data_TemporalSubSpace.trial{1}, sum(rank_source), sobi_delay) ;
            TICA.trial{1} = trial ;
            TICA.unmixing = pinv(Mixing) ; 
    else
        cfg2.method = cfg.method ;
        [TICA] = ft_componentanalysis(cfg2, data_TemporalSubSpace) ;                    % Running ICA for the temporal subsapce
    end
else
        [TICA] = ft_componentanalysis(cfg2, data_TemporalSubSpace) ;                    % Running ICA for the temporal subsapc
end
    
   

%G = U_D*comp.unmixing;
Mixing = pinv(TICA.unmixing);                                   
SpatialICs = SpatialPCs*Mixing ;                                                % Calculating the spatial maps of the ICs

for I = 1:sum(rank_source)
    TICA.label{I} = strcat('IC',num2str(I));
end
%% Restructuring the data

TemporalPCs_FsOrig = data_TemporalSubSpace ; 
TemporalPCs_FsOrig = rmfield(TemporalPCs_FsOrig, 'Freqs') ;
TemporalPCs_FsOrig.continuous = TemporalPCs_FsOrig.trial ;


TICA_FsOrig = TICA;
TICA_FsOrig.continuous{1} = TICA_FsOrig.trial{1} ; 
rmfield(TICA_FsOrig, 'trial') ;
for Trial_Index = 1:length(data_short.trial)
    TICA_FsOrig.trial{Trial_Index}        = TICA_FsOrig.continuous{1}(:,(Trial_Index - 1)*size(data_short.trial{Trial_Index},2) + 1        : Trial_Index*size(data_short.trial{Trial_Index},2)) ; 
    TemporalPCs_FsOrig.trial{Trial_Index} = TemporalPCs_FsOrig.continuous{1}(:,(Trial_Index - 1)*size(data_short.trial{Trial_Index},2) + 1 : Trial_Index*size(data_short.trial{Trial_Index},2)) ; 
end

if  isfield(data_short, 'time') 
    TICA_FsOrig.time = data_short.time ; 
    TemporalPCs_FsOrig.time = data_short.time ; 
end
if  isfield(data_short, 'trialinfo') 
    TICA_FsOrig.trialinfo = data_short.trialinfo ; 
    TemporalPCs_FsOrig.trialinfo = data_short.trialinfo ; 
end
if  isfield(data_short, 'sampleinfo') 
    TICA_FsOrig.sampleinfo = data_short.sampleinfo ; 
    TemporalPCs_FsOrig.sampleinfo = data_short.sampleinfo ;
end


%% Producing the scalp maps from the Spatial ICs and PCs
SpatialICs_Maps = zeros(No_Chan, sum(rank_source)) ;
SpatialPCs_Maps = zeros(No_Chan, sum(rank_source)) ;

for Curren_comp = 1:sum(rank_source)
    SpatialICs_Maps(:,Curren_comp) = SpatialICs(1:No_Chan,Curren_comp) ;
    SpatialPCs_Maps(:,Curren_comp) = SpatialPCs(1:No_Chan,Curren_comp) ;
end


%% Returning the processed data
MBICA.TemporalPCs = TemporalPCs_FsOrig ;
MBICA.TemporalICs = TICA_FsOrig ;
MBICA.SpatialICs = SpatialICs ;
MBICA.SpatialPCs = SpatialPCs_Maps ;
MBICA.MixingMatrix = Mixing ;
MBICA.label = dataToRun.label ;



















