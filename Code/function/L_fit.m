function [fitresult, gof] = L_fit(freq, L, outliers)
%CREATEFIT(FREQ,L)
%  Create a fit.
%
%  Data for 'L fit' fit:
%      X Input : freq
%      Y Output: L
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 29-Jun-2019 21:34:45 自动生成


%% Fit: 'L fit'.
[xData, yData] = prepareCurveData( freq, L );

% Set up fittype and options.
ft = fittype( '0*x + L', 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', outliers);
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMaxChange = 1e-8;
opts.DiffMinChange = 1e-13;
opts.Display = 'Off';
opts.Lower = -0.01;
opts.MaxIter = 80000;
opts.Robust = 'Bisquare';
opts.StartPoint = 5e-07;
opts.TolFun = 1e-12;
opts.TolX = 1e-12;
opts.Upper = 0.01;
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'L(initial fit)' );
% h = plot( fitresult, xData, yData, excludedPoints );
% legend( h, 'L vs. freq', 'Excluded L vs. freq', 'L(initial fit)', 'Location', 'NorthEast' );
% % Label axes
% xlabel freq
% ylabel L
% grid on

end
