%% @ Le Xing & @Alex Casson
% Date: 06/12/2023
% This code is for hBET2 EEG analysis
% Pipeline:
% 1. Load raw data
% 2. select data according to usage log ->
% 3. filtering -> 
% 4. artifacts detection/removal EEGLAB ->
% 5. calculate PSD


%% Intialise Matlab
clear
close all
clc
py.importlib.import_module('fooof'); % assumes Python link is set up correctly


%% Set data
addpath('./libraries/eeglab_current/eeglab2023.1');
addpath('./libraries/fooof_mat-main/fooof_mat-main/fooof_mat');
addpath('./data/hBET 2/');
addpath('./data/hBET 3/');
addpath('./data/hBET 4/');
addpath('./data/hBET 5/');
addpath('./data/hBET 7/');
addpath('./data/hBET 9/');
addpath('./data/hBET 11/');
addpath('./data/hBET 12/');
addpath('./data/hBET 15/');
addpath('./data/hBET 16/');
addpath('./data/hBET 17/');
addpath('./data/hBET 19/');
start_and_end_times = './data/Start and end times.xlsx';


%% Settings
% Note there are some hard coded ASR settings in the pop_clean_rawdata function in calculate_psd
% pop_rejcont settings are also hard coded in calculate_psd

settings.analysis_epoch = 2; % seconds for PSD calculations
settings.segment_size = 1; % minute for cutting up duration
settings.saturation_threshold = 499; % uV for detecting when channel is saturated
settings.artifact_method = 'asr'; % asr or pop_rejcont, will default to asr
settings.k = 10; % ASR setting
settings.wl = 2; % ASR setting
settings.fu = 52.5; % notch filter cut-offs in Hz
settings.fd = 47.5; % notch filter cut-offs in Hz
settings.fl = 50;   % low pass filter cut-off in Hz
settings.fh = 0.16; % higher pass filter cut-off in Hz
settings.n_fft = 2^15; % N for PSD calculations
settings.overlap = 0; % window overlap for PSD calcualtions
settings.alpha_range_start = 8; % Hz
settings.alpha_range_end = 13;  % Hz
settings.fooof_alpha_only = 'on'; % reject fooof fit based on the alpha band only. on or off. Will default to on
settings.fooof_f_range = [1, 30]; % Hz
settings.fooof_r2_threshold = 0.9; % fit value to threshold against


%% Configurations to consider
%duration = {'limit_to_5min','limit_to_10min','full_duration'}; % sheets in start_and_end_times 
%fooof_usage = {'on','off'}; % whether to use fooof; on or off. Will default to on if option is mis-spelt

duration = {'limit_to_5min'}; % sheets in start_and_end_times 
fooof_usage = {'on'}; % whether to use fooof; on or off. Will default to on if option is mis-spelt



%% Run analysis and produce .mat file results which are then converted to xlsx
for i = 1:length(duration)
    sheet = duration{i};
    for j = 1:length(fooof_usage)
        fooof_flag = fooof_usage{j};
        disp(['Now running: ', sheet, ' with FOOOF: ', fooof_flag]);
        hBET_analysis(start_and_end_times,sheet,settings,fooof_flag);
        disp(' ')
    end
end

