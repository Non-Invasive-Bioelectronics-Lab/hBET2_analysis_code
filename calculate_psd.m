function [psd_integral_alpha, artifact_ratio, r2] = calculate_psd(eeg_during_experiment,f_samp,settings,fooof_usage)


    %% Filter data
    % Note Dreem has in-built filter from 0.4 to 35 Hz
    
    % Notch filter
    fu = settings.fu; fd = settings.fd; % cut-offs in Hz
    wu = fu / (f_samp/2);
    wd = fd / (f_samp/2);
    [bn,an]=butter(1,[wd wu], 'stop');
    eeg_during_experiment = filtfilt(bn,an,eeg_during_experiment);
    
    % Low pass filter
    fl = settings.fl; % Hz
    wl = fl / (f_samp/2);
    [bl,al]=butter(1,wl);
    eeg_during_experiment = filtfilt(bl,al,eeg_during_experiment);
    
    % High pass filter
    fh = settings.fh; % Hz
    wh = fh / (f_samp/2);
    [bh,ah]=butter(1,wh,'high');
    eeg_during_experiment = filtfilt(bh,ah,eeg_during_experiment);


    %% Load data into EEGLAB
    
    % eeglab window will pop up, no actions needed, please ignore
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
    EEG = pop_importdata('data', eeg_during_experiment', 'dataformat','array', 'xmin',0, 'srate',f_samp, 'nbchan', size(eeg_during_experiment, 2));

    
    %% Use EEGLAB to remove artefacts
    if strcmp(settings.artifact_method,'pop_rejcont')       
        % Method 1: pop_rejcont()  
        analysis_epoch = settings.analysis_epoch;
        [corrected_eeg, selected_regions] = pop_rejcont(EEG,'freqlimit',[35 128] ,'threshold',10, 'epochlength',analysis_epoch, 'overlap',0.25, 'correct','remove', 'contiguous',1, 'addlength',0, 'taper','hamming', 'verbose','off');
        eeg_no_artefacts = corrected_eeg.data';

    else
        % Method 2: ASR artifact removal--- EEGLAB command line
        corrected_eeg = pop_clean_rawdata(EEG, 'FlatlineCriterion', 'off', 'ChannelCriterion', 'off',...
                                            'LineNoiseCriterion', 'off', 'Highpass', 'off',...
                                            'BurstCriterion', settings.k, 'WindowCriterion', settings.wl, ...
                                            'BurstRejection', 'on', 'Distance', 'Euclidian',...
                                            'WindowCriterionTolerances', [-Inf 7] ); 
        eeg_no_artefacts = corrected_eeg.data';
    end


    %% calculate the proportion of artifacts:
    artifact_ratio = 1 - length(eeg_no_artefacts)/length(eeg_during_experiment);


    %% Calculate PSD for the whole recording
    epoch_duration = settings.analysis_epoch; % seconds, lowest frequency to see is 5 Hz, so 0.2 s period, so 2s to see 10 cycles
    n_fft = settings.n_fft; % targetting 0.01 Hz freqeuncy steps, frequency space = f_samp / N, so N > 25000 for f_samp 250, 2^15 is the closest power of 2
    windows= f_samp * epoch_duration;
    overlap = settings.overlap;
    [pxx,f] = pwelch(eeg_no_artefacts, windows, overlap, n_fft, f_samp);
    
    % figure;
    % semilogx(f,pxx,'LineWidth',2);
    % % plot(f,pxx,'LineWidth',2)
    % axis([0.5 100 ylim]); xlabel('Frequency / Hz'); ylabel('Power Spectral Density / dB/Hz')
    
    % Extract alpha power
    % specify the alpha wave range you want to check
    alpha_range_start = settings.alpha_range_start;   % e.g. from 8~Hz
    alpha_range_end = settings.alpha_range_end;    % e.g. to 13~Hz
    idx = find(f>alpha_range_start & f<alpha_range_end); % extract the values inside the alpha range 8-13~Hz
    uncorrected_psd_integral_alpha = 2* trapz(linspace(alpha_range_start, alpha_range_end, length(idx)),pxx(idx(1):idx(end),:)); % calculate the trapz integration in alpha range. x2 to account for energy in the negative frequencies
    uncorrected_psd_integral_alpha = 10*log10(uncorrected_psd_integral_alpha);


    %% FOOOF
    if strcmp(fooof_usage,'off')
        r2 = NaN;
        % do nothing
    else
    
        % Transpose, to make inputs row vectors
        fooof_freqs = f';
        channels = size(eeg_during_experiment, 2);
        for k = 1:channels
            fooof_psd = pxx(:,k)';
        
            % FOOOF settings
            fooof_settings = struct();  % Use defaults
            f_range = settings.fooof_f_range;
        
            % Run FOOOF
            fooof_results = fooof(fooof_freqs, fooof_psd, f_range, fooof_settings, true);
        
            % Print out the FOOOF Results
            %fooof_results
            %fooof_plot(fooof_results)
    
            % Work out r2 fit value in the alpha band
            idx = ( fooof_results.freqs > alpha_range_start ) & (fooof_results.freqs < alpha_range_end);
            cc = corrcoef(fooof_results.power_spectrum(idx),fooof_results.fooofed_spectrum(idx));
            r2_alpha(k) = cc(1,2);
    
            fooof_all_channels(:,k) = (fooof_results.fooofed_spectrum - fooof_results.ap_fit)';
            fooof_all_f = fooof_results.freqs; % assumes all have the same output frequencies
        end
        
        % Do integration for alpha power
        idx = find(fooof_all_f>alpha_range_start & fooof_all_f<alpha_range_end); % extract the values inside the alpha range 8-13~Hz
        fooof_psd_integral_alpha = 2* trapz(linspace(alpha_range_start, alpha_range_end, length(idx)),fooof_all_channels(idx(1):idx(end),:)); % calculate the trapz integration in alpha range. x2 to account for energy in the negative frequencies
        fooof_psd_integral_alpha = 10*log10(fooof_psd_integral_alpha);
    
        % Record r2 values
        if strcmp(settings.fooof_alpha_only,'off')
            r2 = fooof_results.r_squared; % whole spectrum
        else
            r2 = r2_alpha; % alpha only
        end
        r2_threshold = settings.fooof_r2_threshold;
        reject = r2 < r2_threshold;
        fooof_psd_integral_alpha(reject) = NaN;
    end

    %% Select wanted output
    if strcmp(fooof_usage,'off')
        psd_integral_alpha = uncorrected_psd_integral_alpha;
    else
        psd_integral_alpha = fooof_psd_integral_alpha;
    end

