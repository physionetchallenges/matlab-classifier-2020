function [score, label] = get_12ECG_score(data,header_data,classes, model)

    num_classes = length(classes);

    label = zeros([1,num_classes]);
    score = ones([1,num_classes]);
    
    % Use your classifier here to obtain a label and score for each class.
    features = get_12ECG_features(data,header_data);
    [tmp_label,tmp_score] = predict(model,features);

    label(tmp_label)=1;
    score=tmp_score;

end



