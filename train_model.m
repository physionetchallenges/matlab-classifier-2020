function train_model(input_directory, output_directory)

if ~exist(output_directory, 'dir')
    mkdir(output_directory)
end

disp('Running training code...')
train_12ECG_classifier(input_directory,output_directory);

disp('Done.')
end