function driver(input_directory, output_directory)

    
	% Find files.
    files = {};
    for f = dir(input_directory)'
        if exist(fullfile(input_directory, f.name), 'file') == 2 && f.name(1) ~= '.' && all(f.name(end - 2 : end) == 'mat')
            files{end + 1} = f.name;
        end
    end

    if ~exist(output_directory, 'dir')
        mkdir(output_directory)
    end

    % Load model.
    disp('Loading 12ECG model...')
    model = load_12ECG_model();

    % Iterate over files.
    disp('Predicting 12ECG labels...')
    num_files = length(files);
    for i = 1:num_files
        disp(['    ', num2str(i), '/', num2str(num_files), '...'])

        % Load data.
        file_tmp=strsplit(files{i},'.');
        input_file = fullfile(input_directory, file_tmp{1});
        [data,hea_data] = load_challenge_data(input_file);
        [classes,current_score,current_label] = get_12ECG_score(data,hea_data,model);

        save_challenge_predictions(output_directory,files{i}, current_score, current_label,classes);
	
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

function save_challenge_predictions(output_directory,file, scores, labels,classes)

	filename = strsplit(file,'.');
	output_file = ([output_directory filesep filename{1} '.csv']);

	commaHeader = [classes;repmat({','},1,numel(classes))]; %insert commaas
	commaHeader = commaHeader(:)';
	textHeader = cell2mat(commaHeader); %cHeader in text with commas
	%write header to file
	fid = fopen(output_file,'w'); 
	fprintf(fid,'#%s\n',filename{1});
	fprintf(fid,'%s\n',textHeader);
	fclose(fid);


	%write data to end of file
	dlmwrite(output_file,labels,'delimiter',',','-append');
        dlmwrite(output_file,scores,'delimiter',',','-append');

end

