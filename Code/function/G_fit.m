function [fitresult, gof] = G_fit(freq, G, outliers)
%CREATEFIT(FREQ,G)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : freq
%      Y Output: G
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 29-Jun-2019 21:15:45 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( freq, G );

% Set up fittype and options.
ft = fittype( 'poly1' );
excludedPoints = excludedata( xData, yData, 'Indices', outliers);
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Lower = [-1e-2 -1e-2];
opts.Robust = 'Bisquare';
opts.Upper = [1e-2 1e-2];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'G(initial fit)' );
% h = plot( fitresult, xData, yData, excludedPoints );
% legend( h, 'G vs. freq', 'Excluded G vs. freq', 'G(initial fit)', 'Location', 'NorthEast' );
% % Label axes
% xlabel freq
% ylabel G
% grid on

end
