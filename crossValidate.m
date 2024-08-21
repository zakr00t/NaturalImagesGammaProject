function Y_ = crossValidate(X, Y, K)
% Returns a K-fold cross-validated predictions matrix Y_, given a design matrix X and target matrix Y
% Currently designed to work with multiple linear regression only, future iterations may introduce a user-defined choice of model
    N = size(X, 1);
    grab = floor(N/K); % Number of test samples isolated per fold
    idxs = 1:N;
    idxs_ = idxs; % Test sample indices, depleted over folds
    Y_ = NaN*Y; % Initializing Y_
    % Fold loop:
    for k = 1:K
        switch k
            case K % On the final fold, the remaining samples (>= grab) become the test samples
                foldIdxs = idxs_;
            otherwise
                foldIdxs = datasample(idxs_, grab, replace=false); % Sampling (grab qty) w/o replacement
        end
        idxs_ = setdiff(idxs_, foldIdxs); % Depletion
        % Training
        Xtrain = X(setdiff(idxs, foldIdxs), :);
        Ytrain = Y(setdiff(idxs, foldIdxs), :);
        C = Xtrain\Ytrain; % Coefficients matrix
        % Test
        Xtest = X(foldIdxs, :);
        Y_(foldIdxs, :) = Xtest*C;
    end
end