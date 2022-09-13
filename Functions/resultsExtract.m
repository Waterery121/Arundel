function [variable] = resultsExtract(Variable_number,Filename1,Filename2)
% Function for extracting results from TELEMAC-2D selafin(*.slf) file.
% This function is a combination of functions getmeshnodes.m and
% extractnodesolution.m which were written by Jon French.
% Input: 
% Variable_number  - range from 1 to 5 corresponding to varibles VELOCITY U, 
%                    VELOCITY V, WATER DEPTH, FREE SURFACE, SCALAR VELOCITY  
% FILENAME1        - path of the result(*.slf) file
% FILENAME2        - path of the file contaning selected points
% Output:
% variable         - the varible field of selected points
% Baichuan Yang, UCL

if nargin == 0
Variable_number = 3;
end
[NodeID,XYZ,XY] = getmeshnodes(Filename1,Filename2);
[sol,m] = extractnodesolution(NodeID,XYZ,Filename1);
variable = sol.NodeResult(:,Variable_number,:);
variable = squeeze(variable);