function  model = train_12ECG_classifier(input_directory,output_directory)

disp('Loading data...')

% Find files.
input_files = {};
for f = dir(input_directory)'
    if exist(fullfile(input_directory, f.name), 'file') == 2 && f.name(1) ~= '.' && all(f.name(end - 2 : end) == 'mat')
        input_files{end + 1} = f.name;
    end
end

% read number of unique classes
classes = get_classes(input_directory,input_files);

num_classes = length(classes);
num_files = length(input_files);
Total_data=cell(1,num_files);
Total_header=cell(1,num_files);


% Iterate over files.
for i = 1:num_files
    disp(['    ', num2str(i), '/', num2str(num_files), '...'])
    
    % Load data.
    file_tmp=strsplit(input_files{i},'.');
    tmp_input_file = fullfile(input_directory, file_tmp{1});
    
    [data,hea_data] = load_challenge_data(tmp_input_file);
    
    Total_data{i}=data;
    Total_header{i}=hea_data;
    
end

disp('Training model..')

label=zeros(num_files,num_classes);

for i = 1:num_files
    
    disp(['    ', num2str(i), '/', num2str(num_files), '...']);
    
    data = Total_data{i};
    header_data = Total_header{i};
    
    tmp_features = get_12ECG_features(data,header_data);
    
    features(i,:)=tmp_features;

    for j = 1 : length(header_data)
        if startsWith(header_data{j},'#Dx')
            tmp = strsplit(header_data{j},': ');
            tmp_c = strsplit(tmp{2},',');
            for k=1:length(tmp_c)
                idx=find(strcmp(classes,tmp_c{k}));
                label(i,idx)=1;
            end
            break
        end
    end
    
    
end

model = mnrfit(features,label,'model','hierarchical');

save_12_ECG_model(model,output_directory,classes);

end

function save_12_ECG_model(model,output_directory,classes)
% Save results.
tmp_file = 'finalized_model.mat';
filename=fullfile(output_directory,tmp_file);
save(filename,'model','classes','-v7.3');


disp('Done.')
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
classes=sort(classes);
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
