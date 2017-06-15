%%%% 05/25/2017
%%%% Input: first dimension of sample_data and sample_label should be
%%%% observations

function [accuracy, sensitivity, specificity] = ...
    lz_lda_leave_two_out_PSF(sample_data, sample_label, DiscrimType)


%%%% initiation
TN = 0; % num_orig0_pred0
FP = 0; % num_orig0_pred1; false alarm; type I error
FN = 0; % num_orig1_pred0; false negative; type II error
TP = 0; % num_orig1_pred1

num_sample = size(sample_data, 1);

parfor iRound = 1: num_sample %#ok<PFUNK>
    for jRound = iRound+1 : num_sample
        
        train_data = sample_data;
        train_labl = sample_label;
        train_data([iRound,jRound], :) = [];
        train_labl([iRound,jRound])    = [];
        
        test_data = sample_data([iRound,jRound], :);
        test_labl = sample_label([iRound,jRound]);
        
        Md1 = fitcdiscr(train_data, train_labl, 'DiscrimType', DiscrimType);
        
        for itest = 1: 2
            [predLabl,score] = predict(Md1, test_data(itest,:));
            
            switch char(test_labl(itest))
                case 'Resting-State' % resting-state
                    if strcmp(predLabl, 'Resting-State')
                        TN = TN + 1;
                    elseif strcmp(predLabl, 'Task-Based')
                        FP = FP + 1;
                    end
                case 'Task-Based' % task-based
                    if strcmp(predLabl, 'Resting-State')
                        FN = FN + 1;
                    elseif strcmp(predLabl, 'Task-Based')
                        TP = TP + 1;
                    end
            end
        end
    end
end
accuracy    = (TP + TN) / (TP + TN + FP + FN);
sensitivity = TP / (TP + FN);
specificity = TN / (TN + FP);



