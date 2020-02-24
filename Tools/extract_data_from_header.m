function [recording,Total_time,num_leads,Fs,gain,age_data,sex_data]=extract_data_from_header(header_data);

	tmp_hea = strsplit(header_data{1},' ');
	recording = tmp_hea{1};
	num_leads = str2num(tmp_hea{2});
        Fs = str2num(tmp_hea{3});
	Total_time = str2num(tmp_hea{4})/Fs;
        gain = zeros(1,num_leads);

	for ii=1:num_leads
	        tmp_hea = strsplit(header_data{ii+1},' ');
                tmp_gain=strsplit(tmp_hea{3},'/');
                gain(ii)=str2num(tmp_gain{1});
        end



        for tline = 1:length(header_data)
                if startsWith(header_data{tline},'#Age')
			tmp = strsplit(header_data{tline},': ');
			age_data = str2num(tmp{2});
                elseif startsWith(header_data{tline},'#Sex')
			tmp = strsplit(header_data{tline},': ');
			if strcmp(tmp{2},'Female')
				sex_data = 1;
			else
				sex_data = 0;
			end
		end
	end


end
