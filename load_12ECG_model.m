function model = load_12ECG_model()

        filename='finalized_model.mat';
        A=load(filename);
        model=A.model;

end


