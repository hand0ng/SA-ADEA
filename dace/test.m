%clear;clc;
% X1 = 1:10;
% X2 = X1;
% X3 = X1;
% S = [X1', X2', X3'];
% Y = (X1.^2)'+(X2.^2)'+(X3.^2)';
theta0 = 10*ones(size(S,2),1);
[dmodel, perf] = dacefit(S, Y, 'regpoly2', 'corrgauss', theta0);
[f, df] = regpoly0(S);
y = predictor(S(1,:), dmodel);
