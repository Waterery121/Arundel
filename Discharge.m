function [discharge_s,cumvolume] = Discharge(path_num)
% Function to read water flux derived from BlueKenue and convert to
% cumulative volume by integration along time.
% Input: 
% path_num         - input 1 to read training dataset; input 2 to read test dataset  
% Output:
% discharge_s      - smoothed water flux timeseries of 5 compartments
% cumvolume        - cumulative volume timeseries of 5 compartments
% Baichuan Yang, UCL

if path_num == 1
    p = 'C:\Users\Julian\OneDrive - University College London\Matlab code\SLF\Discharge\TrainCom';
elseif path_num == 2
    p = 'C:\Users\Julian\OneDrive - University College London\Matlab code\SLF\Discharge\TestCom';
end

for i = 1:5
    %Open specified file for reading
    path = [p,num2str(i),'.ts3'];
    fid = fopen(path,'r');
    %skip first 23 header lines
    for j=1:23
        header_line =fgetl(fid);
    end
    %read 2 columns of data until end of file
    j = 0;
    while ~feof(fid)
        j = j+1;
        str = fgetl(fid);   
        s=regexp(str,'\s+');   % find space
        discharge(i,j)=str2num(str(s(end):end));
    end
    fclose(fid);
end


%% smooth
discharge(1,:) = discharge(1,:)*-1;
for i = 1:size(discharge,1)
    discharge_s(i,:) = smooth(discharge(i,:),20);
end
%% clean up noises
if path_num == 2
    for i = 4:5
        d = discharge_s(i,:);
        d(d<0.5) = 0;
        discharge_s(i,:) = d;
    end
    d1 = discharge_s(1,:);
    d1(d1<1) = 0;
    discharge_s(1,:) = d1;
    d3 = discharge_s(3,:);
    d3(d3<1.1)=0;
    discharge_s(3,:)=d3;
end
%% Visualize
t = 0:size(discharge,2)-1;
t = t*2*60;   % convert to second
time = t/3600; % convert to hour
figure()
ti = tiledlayout('flow',"TileSpacing","compact","Padding","compact");

for i=1:5
    nexttile
    plot(time,discharge_s(i,:));
    title(['Compartment ',num2str(i)])
end
ylabel(ti,'Discharge m^3/s');
xlabel(ti,'Time (hour)'); 
%% calculate cumulative volume
cumvolume = cumtrapz(t,discharge_s,2);


figure()
sgtitle('Cumulative volume')
for i=1:5
    subplot(3,2,i)
    plot(t,cumvolume(i,:));
    ylabel('Volume m^3');
    xlabel('Time (second)'); 
end





%% Test
% 
% 
% p = strcat(cd,'\SLF\Discharge\TrainCom3.ts3');
% %Open specified file for reading
% fid = fopen(p,'r');
% %skip first 23 header lines
% for j=1:23
%     header_line =fgetl(fid);
% end
% %read 2 columns of data until end of file
% j = 0;
% while ~feof(fid)
%     j = j+1;
%     str = fgetl(fid);   
%     s=regexp(str,'\s+');   % find space
%     discharge(j)=str2num(str(s(end):end));
% end
% fclose(fid);
% %discharge = discharge*-1;
% ds = smooth(discharge,10);
% 
% t = 0:size(discharge,2)-1;
% t = t*2*60;   % convert to second
% time = t/3600; % convert to hour
% figure()
% 
% plot(time,discharge);
% hold on 
% plot(time,ds);
% hold off
% ylabel('Discharge m^3/s');
% xlabel('Time (h)'); 
% legend('raw','smoothed');
% 
% cumvolume = cumtrapz(t,ds',2);
% figure()
% 
% plot(t,cumvolume);
% ylabel('Volume m^3');
% xlabel('Time (second)'); 










