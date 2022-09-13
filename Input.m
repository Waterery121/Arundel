function [SL] = Input(num,max)
% Function for loading the synthetic surge used as model input and pre-processing. 
% The function is based on telemac_analysis.m written by Jon French.
% Input: 
% num              - to choose the surge file for ANN training and test  
% max              - limit of the timestep of input
% Output:
% SL               - sea level timeseries derived from input file
% Baichuan Yang, UCL

% ANN training input
if num == 1
    FILENAME = strcat(cd,'\Surges\Arun_surge2908.txt');
elseif num == 2
    FILENAME = strcat(cd,'\Surges\1000yrSurge.txt');
end
%Open specified file for reading
fid = fopen(FILENAME,'r');
%skip first 5 header lines
for i=1:5
    header_line =fgetl(fid);
end
%read 2 columns of data until end of file
i = 0;
while ~feof(fid)
    i = i+1;
    T(i) = fscanf(fid, '%g', 1);
    SL(i)= fscanf(fid, '%g\n', 1); % input dataset for ANN
end
fclose(fid);
SL = SL(1:2:end); % resample to 2min interval
SL = SL(1:max); % trim to the desired length

t = 0:length(SL)-1;
t = t*2*60/3600;
figure()
plot(t,SL);
title('Input Boundary Condition');
ylabel('Water Elevation (m)');
xlabel('Time (h)');