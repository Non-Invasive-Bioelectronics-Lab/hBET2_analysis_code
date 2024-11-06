function hBET_analysis(start_and_end_times,sheet,settings,fooof_usage)


    %% load usage log spreedsheet 
    experiment_times = readtable(start_and_end_times, 'Sheet', sheet); 
    total_files = height(experiment_times);
    
    for i=389
    %for i=1:total_files
    
        % Extract filename, start_time, end_time for i-th data
        filename = experiment_times.Recording_filename{i};
        start_minute = experiment_times.Startmin(i); 
        end_minute   = experiment_times.Endmin(i);
    
        if isnan(start_minute) || isnan(end_minute)
            results{i,1}.filename = filename;
            results{i,1}.error = "This dataset is not usable";
            continue;
        end
        
        
        %% load data set
        data = edfread(filename);
        info = edfinfo(filename);
        
        % Extract EEG data
        eeg(:,1) = cell2mat(data.EEGF7_O1);
        eeg(:,2) = cell2mat(data.EEGF8_O2);
        eeg(:,3) = cell2mat(data.EEGF8_F7);
        eeg(:,4) = cell2mat(data.EEGF8_O1);
        eeg(:,5) = cell2mat(data.EEGF7_O2);
        % eeg_copy = eeg;
        
        
        % Set sampling frequency and timebase
        f_samp_array = info.NumSamples/seconds(info.DataRecordDuration);
        f_samp = f_samp_array(2); % assumes EEG is stored in channel 2 in the EEG record
    
    
        %% Select EEG and calculate alpha power for whole selected period
        saturation_threshold = settings.saturation_threshold; % saturation detecion: >499 microVolt
        try
            eeg_during_experiment = eeg((start_minute*60*f_samp)+1:(end_minute*60*f_samp),:);
            [alpha_power_whole, artifact_ratio_whole, r2] = calculate_psd(eeg_during_experiment,f_samp,settings,fooof_usage);
            saturation_ratio = sum(eeg_during_experiment(:,:) >= saturation_threshold) / size(eeg_during_experiment,1); 
        catch
            disp('Error in whole record')
            alpha_power_whole = NaN;
            artifact_ratio_whole = NaN;
            r2 = NaN;
            saturation_ratio = NaN;
        end
    
    
        %% alpha power per minute
        segment_size = settings.segment_size;
        no_of_segments = floor(length(eeg_during_experiment)/(segment_size*f_samp*60));
    
        result_segment = {};
        for j = 1:no_of_segments
    
            try
                eeg_during_segment = eeg_during_experiment( ((j-1)*f_samp*segment_size*60+1 : j*f_samp*segment_size*60), :);
                [alpha_power_segment, artifact_ratio_segment, r2_segment] = calculate_psd(eeg_during_segment,f_samp,settings,fooof_usage);
    
                result_segment{j,1}.alpha_power_segment = alpha_power_segment;
                result_segment{j,1}.artifact_ratio_segment = artifact_ratio_segment;
                if strcmp(fooof_usage,'off')
                    % do nothing
                else
                    result_segment{j,1}.r2 = r2_segment;
                end
    
            catch 
                disp('Error in segment processing')
                result_segment{j,1}.alpha_power_segment = NaN;
                result_segment{j,1}.artifact_ratio_segment = NaN;
                if strcmp(fooof_usage,'off')
                    % do nothing
                else
                    result_segment{j,1}.r2 = NaN;
                end
    
            end
    
        end
    
    
        %% save results
        results{i,1}.filename = filename;
        results{i,1}.alpha_power_whole = alpha_power_whole;
        results{i,1}.artifact_ratio_whole = artifact_ratio_whole;
        results{i,1}.saturation_ratio = saturation_ratio;
        if strcmp(fooof_usage,'off')
            % do nothing
        else
            results{i,1}.r2 = r2;
        end
    
        results{i,1}.no_of_segment = no_of_segments;
        results{i,1}.segment = result_segment;
        results{i,1}.settings = settings;
    
        %% release memory
        clear eeg eeg_during_experiment alpha_power_whole artifact_ratio_whole saturation_ratio r2 eeg_during_segment alpha_power_segment artifact_ratio_segment r2_segment EEG
    
    end
    
    
    %% Save results
    fn = strcat('hBET_results_',sheet,'_fooof_',fooof_usage,'_k',string(settings.k),'_w',string(settings.wl),'.mat');
    save(fn,"results");


