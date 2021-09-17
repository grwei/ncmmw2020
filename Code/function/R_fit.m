function [fitresult, gof] = R_fit(freq, R, outliers)
%CREATEFIT(FREQ,R)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : freq
%      Y Output: R
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 29-Jun-2019 18:19:13 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( freq, R);

% Set up fittype and options.
ft = fittype( 'R0+sqrt(x)*Rs', 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices', outliers);
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMaxChange = 1e-10;
opts.DiffMinChange = 1e-10;
opts.Display = 'Off';
opts.Lower = [0 0];
opts.MaxIter = 40000;
opts.Robust = 'Bisquare';
opts.StartPoint = [0.317099480060861 0.950222048838355];
opts.TolFun = 1e-12;
opts.TolX = 1e-12;
opts.Upper = [2 1];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'R(initial fit)' );
% h = plot( fitresult, xData, yData, excludedPoints );
% legend( h, 'R vs. freq', 'Excluded R vs. freq', 'R(initial fit)', 'Location', 'NorthEast' );
% % Label axes
% xlabel freq
% ylabel R
% grid on

end
