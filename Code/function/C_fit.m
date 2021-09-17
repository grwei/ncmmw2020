function [fitresult, gof] = C_fit(freq, C, outliers)
%CREATEFIT(FREQ,C)
%  Create a fit.
%
%  Data for 'C fit' fit:
%      X Input : freq
%      Y Output: C
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 29-Jun-2019 21:40:45 自动生成


%% Fit: 'C fit'.
[xData, yData] = prepareCurveData( freq, C );

% Set up fittype and options.
ft = fittype( '0*x + C', 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices',outliers);
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMaxChange = 1e-10;
opts.DiffMinChange = 1e-13;
opts.Display = 'Off';
opts.Lower = -0.01;
opts.MaxIter = 40000;
opts.Robust = 'Bisquare';
opts.StartPoint = 0;
opts.TolFun = 1e-12;
opts.TolX = 1e-12;
opts.Upper = 0.01;
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'C(initial fit)' );
% h = plot( fitresult, xData, yData, excludedPoints );
% legend( h, 'C vs. freq', 'Excluded C vs. freq', 'C(initial fit)', 'Location', 'NorthEast' );
% % Label axes
% xlabel freq
% ylabel C
% grid on

end
