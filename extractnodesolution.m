function [sol,m] = extractnodesolution(NodeID,XYZ,FILENAME1)
%EXTRACTNODESOLUTIONS extract solution from specified nodes of Telemac
%results file (Seraphin format)
%
%
% File: matlab_lib_local\Modelling\Telemac\extractnodesolutions.m
% Dependencies:  Telemac Toolbox (telheadr.m)
%
% Author:    Jon French
%
% Revision history
% Version 1.0  2015  JF - created to replace old Fortran programs (e.g.
% telextract etc) used to process Telemac solutions
% Version 1.1  2020  JF - adapted to visualise 1D series
% Version 1.2  2021  JF - adapted for Idealised Deben experiments
% Version 1.3  2021  JF workaround for uigetfile bug
if nargin <= 2
    if ~ispc; menu('Select Telemac geometry/results file','OK'); end
    
    % % UI for Telemac geometry / results file
    [res_file, path]=uigetfile('*.slf;*.sel;*.res','Select Telemac geometry/results file');
    if res_file==0, return, end % user clicked cancel
    FILENAME1 = [path res_file];
end
% % UI for node coordinate file
% [node_file, path]=uigetfile('*.txt','Select Node XY coordinate file');
% if node_file==0, return, end % user clicked cancel
% FILENAME2 = [path node_file];

% FILENAME2 = mesh_node_file;
% 
% read Telemac binary Serafin geometry or solution to extract header info
m = telheadr(FILENAME1);
% 
% % load set of NodeID and X,Y coords
% IDXY = load(FILENAME2);
% NodeID = IDXY(:,1);
% Xn = IDXY(:,2); Yn = IDXY(:,3); %we don't actually use these here

% if nargin == 1
%  figure(figurenum)
%  % plot mesh
%  patch('faces',m.IKLE,'vertices',m.XYZ,'FaceVertexCData',zeros(size(m.XYZ,1),3), ...
%    'FaceColor','none','EdgeColor','k');
%  axis equal
%  axis off
%  hold on
%  % overplot NodeID locations
%  plot(Xn,Yn,'r.','MarkerSize',20);
%  hold off
% end

sol = []; %structure to contain extracted node solutions and other info

%hold results for nodes and variables and steps
NodeResult = zeros(length(NodeID),m.NBV,m.NSTEPS);

% run through timesteps and retain just the data for the nodes we want
% for i=1:5
for i = 1:m.NSTEPS
%  [~,NodeResult(:,:,i)] = telsteprn(m,NodeID,i);
 m = telstepr(m,i);
 NodeResult(:,:,i) = m.RESULT(NodeID,:);
end


%add results for desired nodes
sol.NodeID = NodeID;
% sol.NodeXY = IDXY(:,2:3);
sol.NodeXY = XYZ;
sol.NodeResult = NodeResult;

%add other useful info to sol structure
sol.NBV = m.NBV; %number of variables in solution
sol.title = m.title;
sol.RECV = m.RECV;
sol.NELEM = m.NELEM;
sol.NPOIN = m.NPOIN;
sol.XYZ = m.XYZ;
sol.timestep = m.timestep;
sol.NSTEPS = m.NSTEPS;
sol.DT = m.DT;

FILENAME3 = [FILENAME1(1:length(FILENAME1)-4) '_node_solution.mat'];

save(FILENAME3,'sol','-v7.3');



