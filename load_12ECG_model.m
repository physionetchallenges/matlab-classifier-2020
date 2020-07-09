function model = load_12ECG_model(model_directory)
	
	out_file='finalized_model.mat';
	filename=fullfile(model_directory,out_file);
        A=load(filename);
        model=A;
	
end


