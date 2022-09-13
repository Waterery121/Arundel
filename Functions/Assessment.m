function [R,N] = Assessment(Y,y)
% Function for assessment of simulation against observation
% using RMSE and NSE.
% Input: 
% Y    - observation of several objects
% y    - simulation of several objects
% Output: 
% R    - array of Root Mean Square Error for each object from input
% N    - array of Root Mean Square Error for each object from input
% Baichuan Yang, UCL

for i = 1:size(Y,1)
    % RMSE
    R(i) = sqrt(mean((Y(i,:) - y(i,:)).^2));  
    % NSE
    a = sum((Y(i,:) - y(i,:)).^2);
    b = sum((Y(i,:) - mean(Y(i,:))).^2);
    N(i) = 1 - a/b;
end
Rmean = mean(R);
Nmean = mean(N);
% Display results
fprintf('mean RMSE: %.2f\n',Rmean);
fprintf('mean NSE: %.2f\n',Nmean);