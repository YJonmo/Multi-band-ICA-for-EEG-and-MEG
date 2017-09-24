% This sample code runs the multi-band ICA algorith and the wide-band ICA on the EEG data containing 3 simualtes source of 5 Hz, 12 Hz and 30 HZ with SNRs = 0.1, 0.04, 0.04, respectively.
load('SampleData4MultICA.mat')
cfg_Mul = [] ; 
%% Defining the filter banks
cfg_Mul.bandpass = [ 0.5 4; 3.5 8; 7.5 13; 12 30; 29 100; 47 53 ]       % Defining the filter bands in Hz

%% Defining the number of components for each band
% cfg_Mul.rankThresh = [ 10 ; 10 ; 10; 12; 20;] ;                       % You can use this option if you want to decide the number of components from
% each band as a percentage of the first eigen-value of that band. For example cfg_Mul.rankThresh = [ 10 ; 15 ... ] means that for band one, every
% principal component whose corresponding eigen-value is greater than 10% of the first eigen-value of that band will be in included. This is the same for
% band two but 15% instead of 10%. 

% cfg_Mul.rankThreshTotl = [10] ;                                       % You can use this option if you want to decide the number of components from
% each band as a percentage of the first eigen-value of the band. For example cfg_Mul.rankThreshTotl = [ 5 ] means that for every band, every
% principal component whose corresponding eigen-value is greater than 10% of the first eigen-value of that band will be in includ
cfg_Mul.bandRanks = [ 10; 10; 10; 10; 20; 4] ;                              % Manually define how many principal components from each band should be inlcuded 

%% Duration reduction
% cfg_Mul.trials = 1:20 ;                                               % Definnig how many trial should be sent to the algorithm            
% cfg_Mul.timewindow = [0 45];                                          % Defining the duration of data to be sent to the algorith (only for single
% trila data)

%% Defining the type of PCA algorithm
% If you aim to perform use the algorithm for artifat removal then you neen
% to use cfg_Mul.pca = 'svd'. Otherwise you can leave it  empty and it uses
% eigen decomposition which is faster.
cfg_Mul.pca = 'svd' ;
%% Choosing the ICA algorithm 
% SOBI gives the best decomposition with the multi-ICA algorithm
cfg_Mul.method = 'sobi' ; % Other methods 'runica', 'fastica', 'binica', 'pca', 'svd', 'jader', 'varimax', 'dss', 'cca', 'sobi', 'white'
cfg_Mul.elec.label = SampleData.label ;
[Multi_ICA] = MultiBand_ICA_Simplified(cfg_Mul, SampleData);

%% Now Running the wide-band ICA
cfg_Wide = [] ;
%cfg_Wide.pca = 'svd' ;
cfg_Mul.method = 'sobi' ; % Other methods 'runica', 'fastica', 'binica', 'pca', 'svd', 'jader', 'varimax', 'dss', 'cca', 'sobi', 'white'
cfg_Wide.elec.label = SampleData.label ;
[Wide_ICA] = MultiBand_ICA_Simplified(cfg_Wide, SampleData);




%% Plot the ICA data
cfg = [];
cfg.viewmode = 'vertical';
cfg.channelcolormap = [0 0 0];  %Above two lines plot in blue (events will be black)
cfg.ploteventlabels = 'colorvalue';
cfg.ylim = [-1 1]; %sets yscale
cfg.blocksize = 20;
cfg.continuous              = 'yes';


ft_databrowser(cfg,Multi_ICA.TemporalICs);
ft_databrowser(cfg,Wide_ICA.TemporalICs);
