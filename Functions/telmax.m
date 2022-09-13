function [m, INFILE, OUTFILE] = telmax(INFILE,OUTFILE)
% Example function that reads a Telemac2D (Seraphin) result file
% and outputs the maximum value of the data for all the timesteps to a new 
% Seraphin file
%
% usage: telmax(INFILE,OUTFILE)
%
% INFILE = input filename (optional - a uimenu will appear if no inputs)
% OUTFILE = output filename (optional - an output file will be created if omitted)
%
% other functions that are called:
%
% telheadr.m - for reading Seraphin file header
% telheadw.m - for reading Seraphin file timestep
% telstepr.m - for writing Seraphin file header
% telstepw.m - for writing Seraphin file timestep
%
% modified by JF from telmean.m by Thomas Benson, HR Wallingford
% email: t.benson@hrwallingford.co.uk
% release date: 13-Aug-2009

% check inputs inputs
if nargin<1
    INFILE=[];
end
if nargin<2
    OUTFILE = [];
end

% read the file header
m = telheadr(INFILE);

% get the filename parts
INFILE = m.filename;
[NAME,EXT] = fileparts(INFILE);

% create an averaging array
MAXARRAY = zeros(size(m.RESULT));

fprintf('\nNumber of timesteps: %d\n',m.NSTEPS)
fprintf('Current timestep: %10d',0)

FE = ones(m.NSTEPS,1); %setup flood - ebb sign vector;
%FE = 1 on flood; FE = -1 on ebb
FE(1) = 1; %arbitrary - this will never be the max values so it should 
%not affect the computation of the maxima

% loop through the timesteps
% for i=1:m.NSTEPS
for i=2:m.NSTEPS
    
    fprintf('\b\b\b\b\b\b\b\b\b\b%10d',i);
    
    % read the timestep
    m = telstepr(m,i);

    % add the data to the avergaing array
%     MAXARRAY = max(m.RESULT,MAXARRAY);
    MAXARRAY = ((MAXARRAY+m.RESULT) + abs(MAXARRAY-m.RESULT))/2;

end
fprintf('\n');

% make sure the input file is closed
fclose(m.fid);

% divide by the number of timesteps to get the averages
% and put back into the structure result array
% m.RESULT = AVGARRAY./m.NSTEPS;

m.RESULT = MAXARRAY;

% reset the time to zero for outputting
m.AT = 0;

% if no output file is specified, create a filename
if isempty(OUTFILE)
    OUTFILE = [NAME '_max' EXT];
end

% if the function is called without an output, write to file.
if nargout<1
    % write the header to output
    fid = telheadw(m,OUTFILE);
    % write the timestep (averaged data) to output
    fid = telstepw(m,fid);
    % close the file
    fclose(fid);
end

    