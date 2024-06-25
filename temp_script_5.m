% Info. needed to retrieve datasets:
details.datasets.train.set = "TL";
details.datasets.train.date = '240817';
[details.experimentalDetails, details.datasets.train.idx] = ...
getExperimentalDetails(details.datasets.train.date);

details.datasets.test.set = "AF";
details.datasets.test.date = '010817';
[~, details.datasets.test.idx] = ...
getExperimentalDetails(details.datasets.test.date);

% Hyperparameters:
hparams.radiusMatrixDeg = 0.3:0.3:2;
hparams.patchSizeDeg = 2;
hparams.powerOption = 3; % 1 - EnergyST. 2 - EnergyST/EnergyBL. 3 - Like 2 but with difference between high gamma and gamma ranges.
hparams.feat_combos = ["H", "S", "V", "R", "P", "C", "HS", "HSV", "HSVR", "HSVRPC"]; % The feature combinations for which we want to compute regresion model performance (R^2)
hparams.feat_slices = {1:3, [1, 4], [1, 5], [1, 6], [1, 7], [1, 8], 1:4, 1:5, 1:6, 1:8}; % The columns corresponding to feat_combos. Refer design.mat.
hparams.patch_cutoff = 0.3; % Used for good patch subset selection

% Creating the design tensors for train and test datasets, and obtaining the target matrix (gamma energy in stimulus period):
temp = fieldnames(details.datasets);
for i = 1:length(temp)
    if isfile(temp{i} + ".mat")
        temp_2 = load(temp{i} + ".mat");
        details.datasets.(temp{i}).X = temp_2.X;
        details.datasets.(temp{i}).Y = temp_2.Y;
    else
        [details.datasets.(temp{i}).X, details.datasets.(temp{i}).Y] = ...
        design(details.datasets.(temp{i}).idx, hparams.radiusMatrixDeg,...
        hparams.patchSizeDeg, hparams.powerOption);
    
        P = readmatrix("predict.xlsx", "Sheet", details.datasets.(temp{i}).set);
        P = reshape(P(~isnan(P)), size(P, 1), size(details.datasets.(temp{i}).X, 3));
        P = reshape(P, size(P, 1), 1, size(P, 2)); % Convert to same dimensions as X_train
    
        C = readmatrix("compress.xlsx", "Sheet", details.datasets.(temp{i}).set);
        C = reshape(C(~isnan(C)), size(C, 1), size(details.datasets.(temp{i}).X, 3));
        C = reshape(C, size(C, 1), 1, size(C, 2)); % Convert to same dimensions as X_train
        
        details.datasets.(temp{i}).X = cat(2, details.datasets.(temp{i}).X, P, C);
        temp_2 = details.datasets.(temp{i}); % You simply can't win against MATLAB. You have to make temporary variables. You can't index function calls. Can't save a substruct. Can't enumerate a list. Might as well stop breathing because MATLAB will throw up an error for that too! :)
        save(temp{i} + ".mat", "-struct", "temp_2", "X", "Y")
    end
    temp_2 = cell(1, size(details.datasets.(temp{i}).Y, 2));
    temp_3 = temp_2;
    for j = 1:length(temp_2)
        idx = details.datasets.(temp{i}).X(:, 6, j) > hparams.patch_cutoff;
        temp_2{j} = details.datasets.(temp{i}).X(idx, :, j);
        temp_3{j} = details.datasets.(temp{i}).Y(idx, j);
    end
    % Creating cells to store electrode-wise good patch design matrices and their correseponding target vectors
    details.datasets.(temp{i}).X_sel = temp_2;
    details.datasets.(temp{i}).Y_sel = temp_3;
end
clear i j idx temp_2 temp_3 P C

% Preparing results struct:
for i = 1:length(temp)
    results.(temp{i}).all = [];
    results.(temp{i}).sel = [];
end

% Training linear regression on the training design tensor:
temp = ones(1, length(details.datasets.train.Y_sel));
for idx = 1:length(hparams.feat_combos)
    % All electrodes
    [~, results.train.all.C{idx}, results.train.all.corrs{idx}] = ...
    lr_electrode(details.datasets.train.X(:, hparams.feat_slices{idx}, :), ...
    details.datasets.train.Y);
    results.train.all.corrs_mean(idx) = mean(results.train.all.corrs{idx});
    results.train.all.corrs_SEM(idx) = std(results.train.all.corrs{idx})/sqrt(length(temp));
    % Selected electrodes
    temp_X = details.datasets.train.X_sel;
    for j = 1:length(temp)
        temp_X{j} = details.datasets.train.X_sel{j}(:, hparams.feat_slices{idx});
        temp(j) = length(details.datasets.train.Y_sel{j}) > ...
        length(hparams.feat_slices{idx}); % To prevent trivial fits and overfitting with linear regression. We use up to 7 features at the moment (refer design.mat). An N-1 dimensional/features hyperplane can trivially fit N datapoints. So any less than 9 patches will result in a perfect fit. Can't exceed 32 (total number of samples).
    end
    temp = logical(temp);
    [~, results.train.sel.C{idx}, results.train.sel.corrs{idx}] = ...
    lr_electrode(temp_X, details.datasets.train.Y_sel);
    results.train.sel.C{idx}(:, ~temp) = NaN;
    results.train.sel.corrs{idx}(~temp) = NaN;
    results.train.sel.corrs_mean(idx) = mean(results.train.sel.corrs{idx}, 'omitnan');
    results.train.sel.corrs_SEM(idx) = std(results.train.sel.corrs{idx}, 'omitnan')/sqrt(sum(temp));
end
clear temp temp_X i j idx

% Using coefficients from the trained models for each feature combination to compute corrs on test data:
results.test.all.corrs = results.train.all.corrs; % Initialization
results.test.sel.corrs = results.train.sel.corrs;
for idx = 1:length(hparams.feat_combos)
    results.test.all.Y_{idx} = details.datasets.test.Y; % Initialization
    results.test.sel.Y_{idx} = details.datasets.test.Y_sel;
    for j = 1:size(details.datasets.test.Y, 2)
        % All electrodes
        results.test.all.Y_{idx}(:, j) = ...
        details.datasets.test.X(:, hparams.feat_slices{idx}, j)*...
        results.train.all.C{idx}(:, j);
        temp = corrcoef(details.datasets.test.Y(:, j), results.test.all.Y_{idx}(:, j));
        results.test.all.corrs{idx}(j) = temp(1, 2);
        % Selected electrodes
        results.test.sel.Y_{idx}{j} = ...
        details.datasets.test.X_sel{j}(:, hparams.feat_slices{idx})*...
        results.train.sel.C{idx}(:, j);
        temp = corrcoef(details.datasets.test.Y_sel{j}, results.test.sel.Y_{idx}{j});
        results.test.sel.corrs{idx}(j) = temp(1, 2);
    end
    % All electrodes
    results.test.all.corrs_mean(idx) = mean(results.test.all.corrs{idx});
    results.test.all.corrs_SEM(idx) = std(results.test.all.corrs{idx})/sqrt(size(details.datasets.test.Y, 2));
    % Selected electrodes
    results.test.sel.corrs_mean(idx) = mean(results.test.sel.corrs{idx}, 'omitnan');
    results.test.sel.corrs_SEM(idx) = std(results.test.sel.corrs{idx}, 'omitnan')/sqrt(sum(~isnan(results.test.sel.corrs{idx})));
end
clear temp

% Modification of runAnalyzeData to obtain the tuning model equivalents of the above:
temp = fieldnames(results);
posList = [details.datasets.train.idx, details.datasets.test.idx];
for i=1:length(posList)
    tmp = details.experimentalDetails{posList(i)};
    
    subjectName = tmp{1};
    imageFolderName = tmp{4};
    expDate = tmp{5};
    protocolName = tmp{6};
    if tmp{2} == "Humans"
        imageIndices = 1:16;
    else
        imageIndices = 1:32;
    end

    disp(['Working on ' subjectName expDate protocolName]) % ', set: ' dataType]);
    [results.(temp{i}).all.tuning_corrs, results.(temp{i}).sel.tuning_corrs, ~, predictionString] = ...
        analyzeData(subjectName, expDate, protocolName, imageFolderName, ...
        imageIndices, hparams.powerOption, [], hparams.radiusMatrixDeg);   

    results.(temp{i}).all.tuning_corrs_mean = mean(results.(temp{i}).all.tuning_corrs, 2)';
    results.(temp{i}).all.tuning_corrs_SEM = ...
    (std(results.(temp{i}).all.tuning_corrs, [], 2)/sqrt(size(results.(temp{i}).all.tuning_corrs, 2)))';
    results.(temp{i}).sel.tuning_corrs_mean = mean(results.(temp{i}).sel.tuning_corrs, 2)';
    results.(temp{i}).sel.tuning_corrs_SEM = ...
    (std(results.(temp{i}).sel.tuning_corrs, [], 2)/sqrt(size(results.(temp{i}).sel.tuning_corrs, 2)))';

end
clear tmp subjectName expDate imageFolderName protocolName posList imageIndices ...
    experimentalDetails % Too many annoying variables!!

% We need to compare regression and tuning side-by-side:
predictionString = string(predictionString); 
idx = logical(ones(size(hparams.feat_combos)));
tmp = setdiff(hparams.feat_combos, predictionString);
for k = 1:length(tmp)
    if ~ismember(tmp(k), predictionString)
        idx(hparams.feat_combos == tmp(k)) = 0;
    end
end
clear tmp predictionString

temp_2 = fieldnames(results.train);
temp_3 = ["tuning_corrs_mean", "tuning_corrs_SEM"];
for i = 1:length(temp)
    for j = 1:length(temp_2)
        for k = 1:length(temp_3)
            temp_4 = NaN*zeros(size(hparams.feat_combos));
            temp_4(idx) = results.(temp{i}).(temp_2{j}).(temp_3{k});
            results.(temp{i}).(temp_2{j}).(temp_3{k}) = temp_4;
        end
    end
end
clear temp temp_2 temp_3 temp_4 idx i j k

% Seeing our results so far:
figure 
sgtitle('Comparison of Regression and Tuning Model Performance')
rows = fieldnames(results);
cols = fieldnames(results.train);
for i = 1:length(rows)
    for j = 1:length(cols)
        subplot(2, 2, ((i - 1)*length(cols) + j))
        bar_mean = [results.(rows{i}).(cols{j}).corrs_mean; ...
            results.(rows{i}).(cols{j}).tuning_corrs_mean]';
        bar_SEM = [results.(rows{i}).(cols{j}).corrs_SEM;
            results.(rows{i}).(cols{j}).tuning_corrs_SEM]';
        hBar = bar(bar_mean); % Bar plot handle
        legend('Regression', 'Tuning', 'Location', 'northwest')
        ctr = zeros(size(bar_mean));
        for k = 1:size(bar_mean, 2) % Column index corresponds each bar plot in hBar
            ctr(:, k) = hBar(k).XData + hBar(k).XOffset;   % Note:  XOffset  Is An Undocumented Feature, This Selects The +bar  Centres
        end
        hold on
        errorbar(ctr, bar_mean, bar_SEM, 'ko', 'HandleVisibility','off')
        title(rows{i} + " " + cols{j})        
        xticklabels(hparams.feat_combos)
        xlabel("Feature Combination")
        ylim([0, 1]);
        ylabel("rho")
    end
end
clearvars -except details hparams results