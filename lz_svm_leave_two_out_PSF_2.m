%%%% 05/25/2017
%%%% Input: first dimension of sample_data and sample_label should be
%%%% observations
%%%% modified using function "svmtrain" and "svmclassify" to accomodate 3d
%%%% plots

function [accuracy, sensitivity, specificity] = ...
    lz_svm_leave_two_out_PSF_2(sample_data, sample_label, kernel, boxconstraint)

%%%% initiation
TN = 0; % num_orig0_pred0
FP = 0; % num_orig0_pred1; false alarm; type I error
FN = 0; % num_orig1_pred0; false negative; type II error
TP = 0; % num_orig1_pred1

num_sample = size(sample_data, 1);

parfor iRound = 1: num_sample
    for jRound = iRound+1 : num_sample
        
        train_data = sample_data;
        train_labl = sample_label;
        train_data([iRound,jRound], :) = [];
        train_labl([iRound,jRound])    = [];
        
        test_data = sample_data([iRound,jRound], :);
        test_labl = sample_label([iRound,jRound]);
        
        if strcmp(kernel, '@lz_svm_sigmoid')
            Md1 = svmtrain(train_data, train_labl, 'kernel_function', ...
                @lz_svm_sigmoid, 'autoscale', false, ...
                'boxconstraint', boxconstraint);
        elseif strcmp(kernel, 'polynomial-4')
            Md1 = svmtrain(train_data, train_labl, 'kernel_function', ...
                'polynomial', 'polyorder', 4, 'autoscale', false, ...
                'boxconstraint', boxconstraint);
        elseif strcmp(kernel, 'rbf')
            Md1 = svmtrain(train_data, train_labl, 'kernel_function', ...
                'rbf', 'rbf_sigma', .85, 'autoscale', false, ...
                'boxconstraint', boxconstraint);
        else
            Md1 = svmtrain(train_data, train_labl, 'kernel_function', ...
                kernel, 'autoscale', false, 'boxconstraint', boxconstraint);
        end
        
        for itest = 1: 2
            predLabl = svmclassify(Md1, test_data(itest,:));
            
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


