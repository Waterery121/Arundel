% The main script for validation of TELEMAC-2D model and 
% assessment of ANN emulation model for the Arundel flood 
% inundation model project. The codes can generate several 
% figures for display in the dissertation. Some codes were 
% derived and modified from AnalyseTelemacFloodExtent.m 
% written by Jon French, UCL. And function telmax.m is used 
% which is modified by Jon French as well.
% Dependencies:
% Assessment.m      - to assess the simulation by RMSE and NSE
% Input.m,          - to load input data (tidal timeseries)
% telmax.m          - to read a Telemac2D result file and output 
%                     the maximum value of the data for all the timesteps 
% Discharge.m       - to load discharge of 5 compartments 
%                     calculate and load cumulative volume
% resultsExtract.m  - to load results of given points from 
%                   - TELEMAC-2D model
% Catution:
% The code includes computation of large 3d matrix, 
% RAM of the device should no less than 32G.
% Baichuan Yang, UCL

tic

%%
clc
clear
close all

%% ------ load data ------ %%
% load input and output data
[a,Y] = Discharge(1);     % output of training
[b,Y1] = Discharge(2);    % output of test
P = Input(1,size(Y,2));   % input of training
P1 = Input(2,size(Y1,2)); % input of test
%%
% load trained networks
c1 = load('01_SavedModel\Com1Model0.91.mat').model;
c2 = load('01_SavedModel\Com2Model0.91.mat').model;
c3 = load('01_SavedModel\Com3Model0.95.mat').model;
c4 = load('01_SavedModel\Com4Model0.98.mat').model;
c5 = load('01_SavedModel\Com5Model0.98.mat').model;

Ann = [c1,c2,c3,c4,c5];
for i = 1:5
    closed_net{i} = Ann(i).closed_net; 
    open_net{i} = Ann(i).net; 
end

% simulation
y_open = zeros(size(Y));
y_close = zeros(size(Y));
y_topen = zeros(size(Y1));
y_tclose = zeros(size(Y1));

for i = 1:5
    % Training dataset
    net = open_net{i};
    netc = closed_net{i};
    X = tonndata(P,true,false);
    T = tonndata(Y(i,:),true,false);
    % simulation of training dataset with open-loop
    [x,xi,ai,t] = preparets(net,X,{},T);
    y = net(x,xi,ai);
    y_open(i,(size(Y,2)-size(y,2))+1:end) = cell2mat(y);
    % simulation of training dataset with closed-loop
    [xc,xic,aic,tc] = preparets(netc,X,{},T);
    yc = netc(xc,xic,aic);
    y_close(i,(size(Y,2)-size(yc,2))+1:end) = cell2mat(yc);
    % Test dataset
    X1 = tonndata(P1,true,false);
    T1 = tonndata(Y1(i,:),true,false);
    % simulation of test dataset with open-loop
    [x1,xi1,ai1,t1] = preparets(net,X1,{},T1);
    yt = net(x1,xi1,ai1);
    y_topen(i,(size(Y1,2)-size(yt,2))+1:end) = cell2mat(yt);
    % simulation of test dataset with closed-loop
    [x1,xi1,ai1,t1] = preparets(netc,X1,{},T1);
    ytc = netc(x1,xi1,ai1);
    y_tclose(i,(size(Y1,2)-size(ytc,2))+1:end) = cell2mat(ytc);
end

%% assess the model performance
[RMSE,NSE] = Assessment(Y,y_close);

% Visualize training
figure()
time = 1:size(Y,2);
time = time*2*60/3600;   % convert to hour
tiledlayout('flow','TileSpacing','compact','Padding','compact');
for i=1:5
    nexttile
    plot(time,Y(i,:)/1e5,'k')
    hold on
    txt1 = ['RMSE = ' num2str(round(RMSE(i),0))];
    text('string',txt1,'Units','normalized','position',[0.65,0.15],'FontSize',8);
    txt2 = ['NSE = ' num2str(round(NSE(i),3))];
    text('string',txt2,'Units','normalized','position',[0.65,0.10],'FontSize',8)
    plot(time,y_open(i,:)/1e5,'r')
    plot(time,y_close(i,:)/1e5,'b')
    hold off
    if i==1||i==4;ylabel('Volume (10^5 m^3)');end
    xlabel('Time (hours)')
    title(['Compartment ',num2str(i)])
end
% add title, labels and legend
sgtitle('Training');
legend('boxoff')
legend('TELEMAC-2D','Open-loop NARX','Closed-loop NARX','Fontsize',8,'location','best');

%% Calculate RMSE and NSE
[RMSE,NSE] = Assessment(Y1,y_tclose);

% Test
figure()
time = 1:size(Y1,2);
time = time*2*60/3600;   % convert to hour
tiledlayout('flow','TileSpacing','compact','Padding','compact');
for i=1:5
    nexttile
    plot(time,Y1(i,:)/1e5,'k')
    hold on
    txt1 = ['RMSE = ' num2str(round(RMSE(i),0))];
    text('string',txt1,'Units','normalized','position',[0.65,0.15],'FontSize',8);
    txt2 = ['NSE = ' num2str(round(NSE(i),3))];
    text('string',txt2,'Units','normalized','position',[0.65,0.10],'FontSize',8)
    plot(time,y_topen(i,:)/1e5,'r')
    plot(time,y_tclose(i,:)/1e5,'b')
    hold off
    if i==1||i==4;ylabel('Volume (10^5 m^3)');end
    xlabel('Time (hours)')
    title(['Compartment ',num2str(i)])
end
% add title, labels and legend
sgtitle('Test');
legend('boxoff');
legend('TELEMAC-2D','open-NARX','closed-NARX','Fontsize',8,'location','best');

%% ------ Read DEM and Compartments ------ %%

% Operation on DEM
% read raster file
[Zr,R] = readgeoraster('image/DTM.tif');
info = geotiffinfo('image/DTM.tif');

% assign coordinates to each pixel
Xr = zeros(size(Zr));
Yr = zeros(size(Zr));
for i = 1:size(Zr,1)
    vx = zeros(1,size(Zr,2));
    vy = zeros(1,size(Zr,2));
    for j = 1:size(Zr,2)
        [x,y] = pix2map(info.RefMatrix,i,j);
        vx(j) = x;
        vy(j) = y;
    end
    Xr(i,:) = vx;
    Yr(i,:) = vy;
end

% import outline of the flood plain
filename = strcat(cd,'\shape\outline.shp');
outline = shaperead(filename); % read shape file
lon = extractfield(outline,'X'); % save x coordinates of the vertices
lat = extractfield(outline,'Y'); % save y coordinates of the vertices
% extract Zg by mask
ot_mask = inpolygon(Xr,Yr,lon,lat); 
Zr(~ot_mask) = nan;
%%
% plot DEM
figure()      
surf(Xr/1000,Yr/1000,Zr)
view(2)
%title('Floodplain')
axis equal
shading interp
caxis([0 20])
c = colorbar;
c.Label.String = 'Elevation (m)';

% extract each compartment
outlinefile = strcat(cd,'\shape\_om');
for i = 1:5
    outline = shaperead([outlinefile,num2str(i),'.shp']);
    cx{i} = extractfield(outline,'X');
    cy{i} = extractfield(outline,'Y');
    com_mask = inpolygon(Xr,Yr,cx{i},cy{i});
    Com{i} = Zr;
    Com{i}(~com_mask) = nan;
end
%%
% plot compartments
figure()
t = tiledlayout(2,3,'TileSpacing','Compact','Padding','compact');
for i = 1:5
    nexttile
    surf(Xr/1000,Yr/1000,Com{i})
    view(2)
    title(['Compartment ',num2str(i)])
    shading interp
    caxis([0 3])
    axis equal
end
c=colorbar();
c.Layout.Tile = 'east';
c.Label.String='Elevation (m)';




%% ------ Comparison of Water Depth in controlled points ------ %%

% calculate water depth

D_sim = {};

for i = 1:4
    volume = y_tclose(i,:); % volume of each compartment
    volume(volume<0) = 0;
    % total area, cell size 2*2
    amask = Com{i};
    amask(~isnan(amask)) = 1;
    amask(isnan(amask)) = 0;
    area = sum(sum(amask))*2*2;
    % remove unrelated areas
    Com_tr = Com{i}(~all(isnan(Com{i}),2),:);
    Com_tr = Com_tr(:,~all(isnan(Com_tr),1));
    % average elevation
    Com1d = reshape(Com_tr, size(Com_tr,1)*size(Com_tr,2),1);
    Hmean = mean(Com1d,'omitnan');
    
    % create 3d arrays  
    D = zeros(size(Com_tr,1),size(Com_tr,2),length(volume));
    volume3d = zeros(size(D));
    % The 3rd dimension is variation of volume
    for j=1:length(volume)
        volume3d(:,:,j) = volume(j);
    end
    % calculate water depth in compartment
    D = volume3d/area + Hmean - Com_tr;
    % V<=0 D=0
    D(D<0) = 0;
    D(:,:,volume<=0) = 0;
    D_sim{i} = D;
end

% Com5 is too big for next step 3d arrays 
% Thus, split Com5 to 2 parts vertically
% volume of Compartment5
volume = y_tclose(5,:); 
volume(volume<0) = 0;   
% remove unrelated areas
Com_tr = Com{5}(~all(isnan(Com{5}),2),:);
Com_tr = Com_tr(:,~all(isnan(Com_tr),1));
% split to 2 parts 
Com5{1} = Com_tr(1:720,:);
Com5{2} = Com_tr(721:end,:);
% compute area of each part
amask = Com{5};
amask(~isnan(amask)) = 1;
amask(isnan(amask)) = 0;
area = sum(sum(amask))*2*2; % total area
% calculate mean elevation
Com1d = reshape(Com{5}, size(Com{5},1)*size(Com{5},2),1);
Hmean = mean(Com1d,'omitnan'); % average elevation

for i = 1:2   
    % create 3d arrays  
    D = zeros(size(Com5{i},1),size(Com5{i},2),length(volume));
    volume3d = zeros(size(D));
    for j=1:length(volume)
        volume3d(:,:,j) = volume(j);
    end
    % calculate water depth in compartment
    D = volume3d/area + Hmean - Com5{i};
    D(D<0) = 0;    
    D(:,:,volume==0) = 0;
    D_sim{i+4} = D;
end



%% extract simulations from Telemac of controlled points
FILENAME1 = strcat(cd,'\SLF\Validation.slf');
FILENAME2 = strcat(cd,'\Points\Controled points.txt');
WD = resultsExtract(3,FILENAME1,FILENAME2); % Telemac results of selected points

% indices of the leftright corner of Coms (Oy,Ox)
for i = 1:5
    Oy(i) = find(~all(isnan(Com{i}),2),1);
    Ox(i) = find(~all(isnan(Com{i}),1),1);
end
cxy = load(FILENAME2); % x,y coordinates of controlled locations
Xp = cxy(:,1); Yp = cxy(:,2);
Cind = zeros(length(Xp),3);
% find the nearst cells relative to controlled points
for i = 1:length(Xp)
    xdist = abs(Xp(i) - Xr(1,:)'); 
    xind = find(xdist == min(xdist)); 
    ydist = abs(Yp(i) - Yr(:,1));
    yind = find(ydist == min(ydist)); 
    % the index of nearest cell is (yind,xind)

    % find which comparment the points locate in
    for j = 1:5
        inCom(j) = inpolygon(Xp(i),Yp(i),cx{j},cy{j}); 
    end
    iC = find(inCom==1); % index of the located compartment 
    % indices of controlled points in compartments
    Cind(i,:) = [yind-Oy(iC)+1,xind-Ox(iC)+1,iC];
    if iC == 5
        if Cind(i,1) > 720 % locate in Com5 part2
            Cind(i,3) = 6;
            Cind(i,1) = Cind(i,1)-720;
        end
    end
    D = D_sim{Cind(i,3)};
    d = D(Cind(i,1),Cind(i,2),:);
    % simulated water depth of controlled points
    wd(i,:) = reshape(d,[1,length(volume)]);
end

% Calculate RMSE and NSE
for i = 1:size(wd,1)
    wd(i,1:618) = 0;
end
%%
% assess the model performance
[RMSE,NSE] = Assessment(WD,wd);
%
% Visualize
figure()
time = 1:size(WD,2);
time = time*2*60/3600;
ti = tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
for i=1:size(WD,1)
    nexttile
    plot(time,WD(i,:))
    txt1 = ['RMSE = ' num2str(round(RMSE(i),3))];
    text('string',txt1,'Units','normalized','position',[0.7,0.2],'FontSize',8);
    txt2 = ['NSE = ' num2str(round(NSE(i),3))];
    text('string',txt2,'Units','normalized','position',[0.7,0.1],'FontSize',8)
    hold on
    plot(time,wd(i,:))
    hold off
    if i==4
        legend('boxoff');
        legend('TELEMAC-2D','ANN','Fontsize',8,'location','best');
    end
    
    if Cind(i,3)<6
        title(['Point ',num2str(i),' in Compartment ',num2str(Cind(i,3))])
    else
        title(['Point ',num2str(i),' in Compartment 5'])
    end
end
ylabel(ti,'Water Depth (m)');
xlabel(ti,'Time (hours)');
% Add title, labels and legend
%sgtitle('Comparison of water depth')



%% ------ Comparison of Maximum Flood Extent ------ %%

% Modified from AnalyseTelemacFloodExtent.m
% which was coded by Jon French modified from AnalyseTelemacFPE (Feb 2022)
% Parameters and run control options
VAR = 3;                    %(m) water depth in Selafin file
depth_tol = 0.1;           %(m) used to determine limit of inundation
deepwater_threshold = 2.5;  %(m) used to mask out deep estuary and dock area
%% Call TELMAX to produce a flat file of maximum values
[m, INFILE, OUTFILE] = telmax(FILENAME1,[]);
%extract required data variable
data = m.RESULT(:,VAR);

%% Grid the mesh data for easier processing
XY = m.XYZ;
X = XY(:,1);
Y = XY(:,2);
Z = data;

% interpolate Znode onto grid
Zg = griddata(X,Y,Z,Xr,Yr,'cubic');
Zg(~ot_mask) = nan;
Zg(Zg > deepwater_threshold) = 0;
%% Obtain max flood extent from ANN simulation
Floodmap = zeros(size(Zr));
Oy(6) = Oy(5) + 720;
Ox(6) = Ox(5);
for i = 1:6
    Mfe = max(D_sim{i}(:,:,1:795),[],3,'omitnan');
    [Yc,Xc] = find(Mfe);
    Ym = Yc + Oy(i)-1;
    Xm = Xc + Ox(i)-1;
    for j = 1:length(Ym)
        Floodmap(Ym(j),Xm(j)) = Mfe(Yc(j),Xc(j));
    end
end
Floodmap(~ot_mask) = nan;

%% Create some masks for data processing

% river mask
deepwater = Zg > deepwater_threshold;
% dry area (not flooded above tolerance depth)
dry = Zg < depth_tol;
% area outside model domain (coded as NaN in Zg)
external_area = isnan(Zg);

% compute modelled flood extent from TELEMAC, As
mask0 = ones(size(Zg));
mask0(dry) = 0;
mask0(deepwater) = 0;
mask0(external_area) = 0;

As = sum(sum(mask0))*2*2;
As_mask = mask0;
disp(['TELEMAC-2D flood extent, As, = ' num2str(As) ' m^2']);

% compute modelled flood extent from ANN, Ar
deepwater = Floodmap > deepwater_threshold;
dry = Floodmap < depth_tol;
external_area = isnan(Floodmap);
mask1 = ones(size(Floodmap));
mask1(dry) = 0;
mask1(deepwater) = 0;
mask1(external_area) = 0;

Ar = sum(sum(mask1))*2*2;
Ar_mask = mask1;
disp(['ANN flood extent, Ar, = ' num2str(Ar) ' m^2']);

% compute intersection of As and Ar
Asum = As_mask + Ar_mask;
Ao = Asum == 2;
Ao_mask = zeros(size(Zg));
Ao_mask(Ao) = 1;
Ao = sum(sum(Ao_mask))*2*2;
disp(['Overlap of As and Ar, = ' num2str(Ao) ' m^2']);

% compute F statistic (after Yu et al 2016 Environ. Res. Lett. 11 124011)
F = Ao / (As + Ar - Ao);
disp(['F statistic = ' num2str(F) ]);




%% Visualise flood extents
% show maximum flood depth from TELEMAC-2d
figure() 
ti = tiledlayout(3,2,"TileSpacing","compact","Padding","compact");
%set(gcf,'position',[250 300 1200 300])
ax1 = nexttile;
surf(Xr/1000,Yr/1000,Zg)
view(2)
cmp1 = colormap(ax1,parula);
caxis([0 1.5])
c1 = colorbar;
c1.Label.String = 'Max depth (m)';
axis equal
shading interp
title('TELEMAC-2D Maximum water depth')
set(ax1,'xlim',[501,505]);

% show maximum flood depth from ANN
ax2 = nexttile;
surf(Xr/1000,Yr/1000,Floodmap)
view(2)
cmp2 = colormap(ax2,'parula');
c2 = colorbar;
caxis([0 1.5])
c2.Label.String = 'Max depth (m)';
axis equal
shading interp
title('ANN Maximum water depth')
set(ax2,'xlim',[501,505]);
% show error map
% import colormap
cm=importdata('./seismic_color_map.txt');
seismic=cm(:,1:3);

ax3 = nexttile;
surf(Xr/1000,Yr/1000,Zg-Floodmap)
view(2)
axis equal
shading interp
caxis([-1.5 1.5])
c3 = colorbar;
cmp3 = colormap(ax3,seismic);
c3.Label.String = 'Error (m)';
title('Sptial distribution of the model error')
set(ax3,'xlim',[501,505]);
linkaxes([ax1 ax2 ax3],'xy')


ax4 = nexttile;
imagesc(Xr(1,:)/1000,Y(:,1)/1000,As_mask)
set(ax4,'dataAspectRatio',[1 6 1])
txt1 = ['Flood extent=' num2str(round(As/1e6,2)) ' km^2'];
text('string',txt1,'Units','normalized','position',[0.5,0.1],'FontSize',8,'Color','w','FontWeight','bold');
title('TELEMAC-2D flood area')
ax5 = nexttile;
imagesc(Xr(1,:)/1000,Y(:,1)/1000,Ar_mask)
set(ax5,'dataAspectRatio',[1 6 1])
txt2 = ['Flood extent=' num2str(round(Ar/1e6,2)) ' km^2'];
text('string',txt2,'Units','normalized','position',[0.5,0.1],'FontSize',8,'Color','w','FontWeight','bold');
title('ANN flood area')
ax6 = nexttile;
imagesc(Xr(1,:)/1000,Y(:,1)/1000,Ao_mask)
set(ax6,'dataAspectRatio',[1 6 1])
txt3 = ['F statistic = ' num2str(round(F,3))];
text('string',txt3,'Units','normalized','position',[0.5,0.1],'FontSize',8,'Color','w','FontWeight','bold');
title('Overlapped flood area')
ax4.YAxis.Visible = 'off';
ax5.YAxis.Visible = 'off';
ax6.YAxis.Visible = 'off';
ax4.XAxis.Visible = 'off';
ax5.XAxis.Visible = 'off';
ax6.XAxis.Visible = 'off';

%%
toc
disp(['Run time: ',num2str(toc/60)],' min');

%% ------ Validation of TELEMAC-2D model ------ %%

% load results
FILENAME1 = strcat(cd,'\TELEMAC\V1_0.01','.slf');
FILENAME2 = strcat(cd,'\TELEMAC\V2_0.015','.slf');
FILENAME3 = strcat(cd,'\TELEMAC\V3_0.02','.slf');
FILENAME4 = strcat(cd,'\TELEMAC\V4_0.025','.slf');
FILENAME5 = strcat(cd,'\TELEMAC\V5_0.03','.slf');
FILENAME6 = strcat(cd,'\TELEMAC\V6_0.035','.slf');
OBS = strcat(cd,'\Points\observation_locations.txt');
friction = 0.01:0.005:0.035;
for i = 1:6
    % Telemac results of selected points
    fs = resultsExtract(4,eval(['FILENAME',num2str(i)]),OBS); 
    FS{i} = fs;
end

% load observations
path = strcat(cd,'/TELEMAC/');
Arundel = load([path,'Arundel-17-06-22-1230start-15min-OD-level.txt']);
Littleham = load([path,'Littlehampton-17-06-22-1230start-15min-OD-level.txt']);
obs = [Arundel,Littleham];

step = min(size(obs,1),size(FS{1},2));
% plot best run
time = 0:15*60:15*60*(step-1); %15min interval convert to second
time = time ./3600; %hours
figure()
t = tiledlayout(6,2,"TileSpacing",'compact');
for i = 1:6
    fs = FS{i};
    Y = obs(1:step,:)';
    y = fs(:,1:step);
    [RMSE,NSE] = Assessment(Y,y);
    
    nexttile
    plot(time,obs(1:step,1))
    hold on
    txt = ['Fr=' num2str(friction(i))];
    text('string',txt,'Units','normalized','position',[0.84,0.25], ...
        'FontSize',6,'FontWeight','bold');
    txt1 = ['RMSE=' num2str(round(RMSE(1),2))];
    text('string',txt1,'Units','normalized','position',[0.84,0.15], ...
        'FontSize',6,'FontWeight','bold');
    txt2 = ['NSE=' num2str(round(NSE(1),2))];
    text('string',txt2,'Units','normalized','position',[0.84,0.05], ...
        'FontSize',6,'FontWeight','bold')
    plot(time,FS{i}(1,1:step))
    hold off
    if i==1
        title('Arundel');
        legend('Obs','Sim','Fontsize',6,'location','best');
        legend('boxoff')
    end
    ylim([-1 3])

    nexttile
    plot(time,obs(1:step,2))
    hold on
    txt = ['Fr=' num2str(friction(i))];
    text('string',txt,'Units','normalized','position',[0.84,0.25], ...
        'FontSize',6,'FontWeight','bold');
    txt1 = ['RMSE=' num2str(round(RMSE(2),2))];
    text('string',txt1,'Units','normalized','position',[0.84,0.15], ...
        'FontSize',6,'FontWeight','bold');
    txt2 = ['NSE=' num2str(round(NSE(2),2))];
    text('string',txt2,'Units','normalized','position',[0.84,0.05], ...
        'FontSize',6,'FontWeight','bold')
    plot(time,FS{i}(2,1:step))
    hold off
    if i==1;title('Littlehampton');end
    ylim([-2 3])
end
ylabel(t,'Water level (m OD)');
xlabel(t,'Time (hours)');















