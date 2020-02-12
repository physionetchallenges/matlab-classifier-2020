function [classes,score, label] = get_12ECG_score(data,hea_data, model)

    classes = {'Normal', 'AF', 'I-AVB', 'LBBB', 'RBBB', 'PAC','PVC','STD','STE'};
    num_of_classes = length(classes);

    label = zeros([1,num_of_classes]);
    score = ones([1,num_of_classes]);
    
    % Use your classifier here to obtain a label and score for each class.
    tmp_features = get_12ECG_features(data,hea_data);
    [tmp_label,tmp_score] = predict(model,tmp_features);

    label(tmp_label)=1;
    score=score*.1;%mp_score;

end



