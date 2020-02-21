function [recording,Total_time,num_leads,Fs,gain]=extract_data_from_header(header_data);

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

	
end
