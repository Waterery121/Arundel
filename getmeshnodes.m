function [NodeID,XYZ,outfile,XY] = getmeshnodes(FILENAME1,FILENAME2)
%GETMESHNODES return mesh node nearest to supplied X,Y coordinates (Telemac)
% 
% Usage:
% [NodeID,XYZ] = getmeshnodes()
% 
% will prompt for:
% 1. Telemac geometry / results file (Seraphin format)
% 2. Desired node coordinate file (text, X and Y coordinates)
% 3. Output file to contain Node IDs and their X and Y coordinates
%
% File: matlab_lib_local\Modelling\Telemac\getmeshnodes.m
% Dependencies:  Telemac Toolbox (telheadr.m)
%   
% Author:    Jon French
%            
% Revision history
% Version 1.0  2013  JF
% Version 1.1  2015  JF - added argument to control plotting option

if nargin == 0
% UI for Telemac geometry / results file
if ~ispc; menu('Select Telemac geometry file','OK'); end
[res_file, path]=uigetfile('*.slf','Select Telemac geometry file');
if res_file==0, return, end % user clicked cancel
FILENAME1 = [path res_file];

% UI for node coordinate file
if ~ispc; menu('Select observation locations file','OK'); end
[node_file, path]=uigetfile('*.txt','Select Node XY coordinate file');
if node_file==0, return, end % user clicked cancel
FILENAME2 = [path node_file];

end

% read Telemac binary Serafin geometry or solution to extract header info
m = telheadr(FILENAME1);
mXYZ = m.XYZ;

% load set of X,Y coords
XY = load(FILENAME2);
Xn = XY(:,1); Yn = XY(:,2);

% plot up mesh
% scrnsz = get(0,'ScreenSize');
% figure('Position',[scrnsz(3)*0.01 scrnsz(4)*0.05 scrnsz(3)*0.98 scrnsz(4)*0.9]);
figure(1)
hold off
patch('faces',m.IKLE,'vertices',m.XYZ,'FaceVertexCData',zeros(size(m.XYZ,1),3), ...
   'FaceColor','none','EdgeColor','k');
axis equal
axis off
hold on

% plot desired coordinates
plot(Xn,Yn,'r.','MarkerSize',20);


% search for nearest node in mesh
NodeID = ones(length(Xn),1);     %mesh node IDs that we want
Nodeoffset = ones(length(Xn),1); %discrepancy between desired and actual mesh XY
tolerance_distance = 20; % tolerance value for reporting of large offsets, m
tolcheck = ones(length(Xn),1);
disp(' ');
disp('GetMeshNodes .... v1.1')
for i = 1: length(Xn)
 xdist = Xn(i) - mXYZ(:,1);
 ydist = Yn(i) - mXYZ(:,2);
 dist = sqrt(xdist.^2+ydist.^2);
 NodeID(i) = find(dist == min(dist));
 Nodeoffset(i) = dist(NodeID(i));
 if Nodeoffset(i) > tolerance_distance
  disp(['*** Warning - Point ' num2str((i)) '; Node ' num2str(NodeID(i)) ' has offset of ' num2str(Nodeoffset(i)) 'm']);
  tolcheck(i) = 0;
 end
end
disp(' ');
valid = tolcheck == 1;
mean_offset = mean(Nodeoffset);
adj_mean_offset = mean(Nodeoffset(valid));
disp(['Mean offset = ' num2str(mean_offset) 'm including nodes exceeding ' num2str(tolerance_distance) 'm tolerance and ' num2str(adj_mean_offset) 'm excluding these nodes']);
disp(' ');

% plot actual nearest node coordinates
plot(mXYZ(NodeID,1),mXYZ(NodeID,2),'b.','MarkerSize',20);
text
[path,name,suffix] = fileparts(FILENAME1);
res_file = [name suffix];
outfile = [res_file(1:length(res_file)-4) '_NodeIDs.txt'];
% write NodeID and actual NodeX, NodeY to file
fid = fopen(outfile,'w');
for i = 1:length(NodeID)
 fprintf('%3.0f %7.0f  %7.1f  %7.1f\n',i, NodeID(i), m.XYZ(NodeID(i),1), m.XYZ(NodeID(i),2)); %echo
 fprintf(fid,'%7.0f  %7.1f  %7.1f\n',NodeID(i), m.XYZ(NodeID(i),1), m.XYZ(NodeID(i),2));
end
fclose(fid);

XYZ = m.XYZ(NodeID,:);
