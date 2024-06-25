function [Y_, C, corrs] = lr_electrode(X, Y)
% OLS Linear Regression for tensor/cell inputs with an electrode-numbered depth axis (might change in the future, currently suits the use case)
    if (class(X) == "double") & (class(Y) == "double")
        C = zeros(size(X, 2), size(Y, 2)); % Coefficients matrix
        Y_ = Y; % Gamma Energy predictions matrix. Just initializing, these values will be overwritten with predictions don't worry.
        corrs = zeros(1, size(Y, 2));
        for j = 1:size(Y, 2)
            C(:, j) = X(:, :, j)\Y(:, j); % Other options include pinv(X)*Y(:, j), linsolve(X, Y(:, j)), fitlm(X, Y(:, j))
            Y_(:, j) = X(:, :, j)*C(:, j);
            tmp = corrcoef(Y(:, j), Y_(:, j));
            corrs(j) = tmp(1, 2);
        end
    elseif (class(X) == "cell") & (class(Y) == "cell")
        C = zeros(size(X{1}, 2), length(Y)); % Coefficients matrix
        Y_ = Y;
        corrs = zeros(1, length(Y));
        for j = 1:length(Y)
            C(:, j) = X{j}\Y{j}; 
            Y_{j} = X{j}*C(:, j);
            tmp = corrcoef(Y{j}, Y_{j});
            corrs(j) = tmp(1, 2);
        end
    end
end