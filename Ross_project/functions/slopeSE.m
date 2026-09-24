function se = slopeSE(x, y)
% slopeSE computes the standard error of the slope from linear regression
%
% Input:
%   x - independent variable (vector)
%   y - dependent variable (vector)
% Output:
%   se - standard error of the slope

    % Ensure column vectors
    x = x(:);
    y = y(:);

    % Remove NaN
    mask = ~(isnan(x) | isnan(y));
    x = x(mask);
    y = y(mask);

    n = numel(x);
    if n < 3
        se = NaN; % not enough points
        return
    end

    % Fit line: slope and intercept
    X = [ones(n,1), x];
    b = X \ y;  % least-squares [intercept; slope]
    yhat = X * b;

    % Residuals
    residuals = y - yhat;
    s2 = sum(residuals.^2) / (n - 2); % variance of residuals

    % Sxx term
    Sxx = sum((x - mean(x)).^2);

    % Standard error of slope
    se = sqrt(s2 / Sxx);
end
