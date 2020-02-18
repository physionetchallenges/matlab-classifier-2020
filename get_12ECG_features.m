function features = get_12ECG_features(data, hea_data)

	% read number of leads, sample frequency and gain from the header.	
	tmp_hea = strsplit(hea_data{1},' ');
	number_of_leads = str2num(tmp_hea{2});
	sample_Fs = str2num(tmp_hea{3});
	gain_lead = zeros([1,number_of_leads]);
	
	for ii=1:number_of_leads
		tmp_hea = strsplit(hea_data{ii+1},' ');
		tmp_gain=strsplit(tmp_hea{3},'/');
		gain_lead(ii)=str2num(tmp_gain{1});
	end

	tmp_label = strsplit(hea_data{16},' ');
	label = tmp_label{1};

%   Loop to get R peaks for all leads but in our sample code
%   We are only using data from lead1
%    for ii=1:number_of_leads:
%        [peaks(ii,:),idx(ii,:)] = detect_peaks(data(ii,:),sample_Fs,gain_lead(ii));
%    end
	
	[peaks,idx]=detect_peaks(data(1,:),sample_Fs,gain_lead(1));

	% mean
	mean_RR = nanmean(idx/sample_Fs*1000);
	mean_Peaks = nanmean(peaks*gain_lead(1));
	% median
	median_RR = nanmedian(idx/sample_Fs*1000);
        median_Peaks = nanmedian(peaks*gain_lead(1));

      	% standard deviation
	std_RR = std(idx/sample_Fs*1000);
	std_Peaks = std(peaks*gain_lead(1));

	%    variance
    	var_RR = var(idx/sample_Fs*1000);
    	var_Peaks = var(peaks*gain_lead(1));

	%   Skewness
    	skew_RR = skewness(idx/sample_Fs*1000);
    	skew_Peaks = skewness(peaks*gain_lead(1));

	%   Kurtosis
    	kurt_RR = kurtosis(idx/sample_Fs*1000);
    	kurt_Peaks = kurtosis(peaks*gain_lead(1));

	features = [mean_RR mean_Peaks median_RR median_Peaks std_RR std_Peaks var_RR var_Peaks skew_RR skew_Peaks kurt_RR kurt_Peaks];
end

function [detected_peaks_values,detected_peaks_indices] = detect_peaks(ecg_measurements,signal_frequency,gain)


%        Method responsible for extracting peaks from loaded ECG measurements data through measurements processing.
%        Michał Sznajder (Jagiellonian University) - technical contact (msznajder@gmail.com)
%        Marta Łukowska (Jagiellonian University)
%        Janko Slavic peak detection algorithm and implementation.
%        https://github.com/c-labpl/qrs_detector
%        https://github.com/jankoslavic/py-tools/tree/master/findpeaks

        filter_lowcut = 0.001;
        filter_highcut = 15.0;
        filter_order = 1;
        integration_window = 30;  % Change proportionally when adjusting frequency (in samples).
        findpeaks_limit = 0.35;
        findpeaks_spacing = 100;  % Change proportionally when adjusting frequency (in samples).
        refractory_period = 240; % Change proportionally when adjusting frequency (in samples).
        qrs_peak_filtering_factor = 0.125;
        noise_peak_filtering_factor = 0.125;
        qrs_noise_diff_weight = 0.25;

        qrs_peak_value = 0.0;
        noise_peak_value = 0.0;
        threshold_value = 0.0;

        % Measurements filtering - 0-15 Hz band pass filter.
        filtered_ecg_measurements = bandpass_filter(ecg_measurements, filter_lowcut, filter_highcut, signal_frequency, filter_order);

	% Derivative - provides QRS slope information.
	differentiated_ecg_measurements=diff(filtered_ecg_measurements);

	% Squaring - intensifies values received in derivative.
        squared_ecg_measurements = differentiated_ecg_measurements.^2;

	% Moving-window integration.
        integrated_ecg_measurements = conv(squared_ecg_measurements, ones([1,integration_window]));

	% Fiducial mark - peak detection on integrated measurements.
        detected_peaks_indices = findpeaks(integrated_ecg_measurements,findpeaks_limit,findpeaks_spacing);

        detected_peaks_values = integrated_ecg_measurements(detected_peaks_indices);

end

function y=bandpass_filter(data, lowcut, highcut, signal_freq, filter_order)
       
%        Method responsible for creating and applying Butterworth filter.
%        :param deque data: raw data
%        :param float lowcut: filter lowcut frequency value
%        :param float highcut: filter highcut frequency value
%        :param int signal_freq: signal frequency in samples per second (Hz)
%        :param int filter_order: filter order
%        :return array: filtered data

	nyquist_freq = 0.5 * signal_freq;
        low = lowcut / nyquist_freq;
        high = highcut / nyquist_freq;
        [b, a] = butter(filter_order, [low high], 'bandpass');
        y = filter(b, a, data);
end


function ind = findpeaks(data,limit,spacing)

%	Janko Slavic peak detection algorithm and implementation.
%        https://github.com/jankoslavic/py-tools/tree/master/findpeaks
%        Finds peaks in `data` which are of `spacing` width and >=`limit`.
%        :param ndarray data: data
%        :param float spacing: minimum spacing to the next peak (should be 1 or more)
%        :param float limit: peaks should have value greater or equal
%        :return array: detected peaks indexes array
	

	if isempty(spacing)
		spacing=1;
	end

        len = length(data);
	x = zeros([1 ,len + 2 * spacing]);
        x(1:spacing) = data(1) - 1.e-6;
        x(end-spacing:end) = data(end-1) - 1.e-6;
	x(spacing:spacing + len-1) = data(1,:);
        peak_candidate = true ([1,len]);

        for s =1:spacing-1
            start = spacing - s;
            h_b = x(start: start + len);  % before
            start = spacing;
            h_c = x(start: start + len);  % central
            start = spacing + s;
            h_a = x(start: start + len);  % after

		for ii = length(peak_candidate)
			if (peak_candidate(ii)==1) && ((h_c(ii) > h_b(ii))) && ((h_c(ii) > h_a(ii)))
				peak_candidate(ii)=true;
			end
		end
	end

        ind = find(peak_candidate);
 	
	if ~isempty(limit)
		ind = ind(find(data(ind)>limit));
	end

end

