% This function calculates RR intervals and SQI values using 
% the HRV toolbox (https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox.git)
% more information can be found: 
% Vest A, Da Poian G, Li Q, Liu C, Nemati S, Shah A, Clifford GD, 
% "An Open Source Benchmarked Toolbox for Cardiovascular Waveform and Interval Analysis", 
%  Physiological measurement 39, no. 10 (2018): 105004. DOI:10.5281/zenodo.1243111; 2018. 

% Then, it finds the fiducial points using ECGkit (https://github.com/marianux/ecg-kit.git)
% More information can be found:
% Demski AJ, Llamedo Soria M. ecg-kit a Matlab Toolbox for Cardiovascular Signal Processing. 
% Journal of Open Research Software. 2016;4(1):e8. DOI: http://doi.org/10.5334/jors.86

% It creates a temporal beat alignment using (https://github.com/Tereshchenkolab/Origin.git)
% and then extract the GEH parameters using (https://github.com/Tereshchenkolab/Global-Electrical-Heterogeneity.git)
% More information can be found:
% Erick Perez-Alday,et al; Importance of the Heart Vector Origin Point Definition for an ECG analysis: 
% The Atherosclerosis Risk in Communities (ARIC) study. Comp Biol Med, Volume 104, January 2019, 
% pages 127-138. https://doi.org/10.1016/j.compbiomed.2018.11.013
% Waks JW, et al. Global Electric Heterogeneity Risk Score for Prediction of Sudden Cardiac Death in the General Population: 
% The Atherosclerosis Risk in Communities (ARIC) and Cardiovascular Health (CHS) Studies. Circulation. 2016;133:2222-2234.
%

function features = get_12ECG_features(data, header_data)

       % addfunction path needed
        addpath(genpath('Tools/'))
        load('HRVparams_12ECG','HRVparams')

	% read number of leads, sample frequency and gain from the header.	

	[recording,Total_time,num_leads,Fs,gain]=extract_data_from_header(header_data);

	HRVparams.Fs=Fs;
        HRVparams.PeakDetect.windows = floor(Total_time-1);
        HRVparams.windowlength = floor(Total_time);

	try

                for i =1:num_leads
                        Lead12wGain(i,:) = data(i,:)* gain(i);
                end


                % median filter to remove bw
                for i=1:num_leads
                        ECG12filt(i,:) = medianfilter(Lead12wGain(i,:)', Fs);
                end

                % convert 12Leads to XYZ leads using Kors transformation
                XYZLeads = Kors_git(ECG12filt);

                VecMag = vecnorm(XYZLeads');


                % Convert ECG waveform in rr intervals
                [t, rr, jqrs_ann, SQIvalue , tSQI] = ConvertRawDataToRRIntervals(VecMag, HRVparams, recording);
                sqi = [tSQI', SQIvalue'];

                % Find fiducial points using ECGKit
                ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
                wavedet_config.setup.wavedet.QRS_detection_only = 0;
                [Fid_pts,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);

                [XYZ_Median,Fid_pts_Median] = Time_coherent_code_github(XYZLeads,Fid_pts,Fs);

                features = GEH_analysis_git(XYZ_Median,Fid_pts_Median,Fs);

	catch
		features = NaN(1,22);
	end

end

