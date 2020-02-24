function driver(input_directory, output_directory)

    
	% Find files.
    input_files = {};
    for f = dir(input_directory)'
        if exist(fullfile(input_directory, f.name), 'file') == 2 && f.name(1) ~= '.' && all(f.name(end - 2 : end) == 'mat')
            input_files{end + 1} = f.name;
        end
    end

    if ~exist(output_directory, 'dir')
        mkdir(output_directory)
    end

    
    % read number of unique classes
    classes = get_classes(input_directory,input_files);


    % Load model.
    disp('Loading 12ECG model...')
    model = load_12ECG_model();

    % Iterate over files.
    disp('Predicting 12ECG labels...')
    num_files = length(input_files);
    for i = 1:num_files
        disp(['    ', num2str(i), '/', num2str(num_files), '...'])

        % Load data.
        file_tmp=strsplit(input_files{i},'.');
        tmp_input_file = fullfile(input_directory, file_tmp{1});
        [data,header_data] = load_challenge_data(tmp_input_file);
        [current_score,current_label] = run_12ECG_classifier(data,header_data,classes,model);

        save_challenge_predictions(output_directory,file_tmp{1}, current_score, current_label,classes);
	
    end


    disp('Done.')
end


function [data,tlines] = load_challenge_data(filename)

        % Opening header file
        fid=fopen([filename '.hea']);
        if (fid<=0)
                disp(['error in opening file ' filename]);
        end

        tline = fgetl(fid);
        tlines = cell(0,1);
        while ischar(tline)
            tlines{end+1,1} = tline;
            tline = fgetl(fid);
        end
        fclose(fid);

        f=load([filename '.mat']);
        try
                data = f.val;
        catch ex
                rethrow(ex);
        end

end

% find unique number of classes
function classes = get_classes(input_directory,files)
	
	classes={};
	num_files = length(files);
	k=1;
    	for i = 1:num_files
		g = strrep(files{i},'.mat','.hea');
		input_file = fullfile(input_directory, g);
	        fid=fopen(input_file);
	        tline = fgetl(fid);
        	tlines = cell(0,1);

		while ischar(tline)
        	    tlines{end+1,1} = tline;
	            tline = fgetl(fid);
			if startsWith(tline,'#Dx')
				tmp = strsplit(tline,': ');
				tmp_c = strsplit(tmp{2},',');
				for j=1:length(tmp_c)
		                	idx2 = find(strcmp(classes,tmp_c{j}));
		                	if isempty(idx2)
                	        		classes{k}=tmp_c{j};
                        			k=k+1;
                			end
				end
			break
        		end
		end
        	fclose(fid);
	end
	classes=sort(classes)
end



function save_challenge_predictions(output_directory,recording, scores, labels,classes)

	output_file = ([output_directory filesep recording '.csv']);

	Total_classes = strjoin(classes,','); %insert commaas
	%write header to file
	fid = fopen(output_file,'w');
	fprintf(fid,'#%s\n',recording);
	fprintf(fid,'%s\n',Total_classes);
	fclose(fid);

	%write data to end of file
	dlmwrite(output_file,labels,'delimiter',',','-append','precision',4);
        dlmwrite(output_file,scores,'delimiter',',','-append','precision',4);

end

