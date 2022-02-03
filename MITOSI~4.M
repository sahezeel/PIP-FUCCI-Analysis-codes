close all 
clear, clc 
tic
%%
%%===================   ===== User inputs ============================
%path to the .csv file:
csv_path = 'C:\Users\sundarer\Desktop\GalbanLab\Sahezeel_Al_MIL\New Test Data for Ram\Red drop Data\10 palblo\';
%subfolder containing graph files, this subfolder should be in the same folder as this matlab file
graph.subdir = 'Mitosis Stuck in G11';
%prefix and suffix of graph file names:
graph.prefix = 'Track ID';
graph.suffix = 'png';

%save names:
SaveGraphName = 'Median_Mitosis_Stuck_in_G1.fig';
%SaveName = 'CrossOverFrames1';

%pass band for lowpass filter:
pass_band = 0.1;%choosing from (0,1)

%xticks
x_ticks = 2;

%window size for moving average filter to smooth the graphs
win = 1;
win2 = 50;
win3 = 100;
% win4 = 250;

%The color of the confidence intervals is decided by setting this value
%between 0 and 1. 0 is blank with no color and 1 is completely rich
patchSatVal = 0.1;

%Low and High values for removing outliers lines from all the lines we have
%in both red and green channel

low = 3.5;
high = 97;

%%
%Preparation:
graph.dir = [csv_path,graph.subdir,'\'];
addpath(graph.dir);
Len_graph_prefix = length(graph.prefix);%length of prefix of graph file names
tracks_file = [csv_path,'tracks.mat'];
load(tracks_file)% tracks;%tracks

%%
%list graph files:
%find max num of frames in all tracks defined by graph files:
[Max_Num_frames,List_of_graph_files] = findMaxNumFrames(graph,tracks);
Num_of_graph_files = length(List_of_graph_files);%number of graph files

%for crossover:
id = [];
id1 = [];
drop_frame = [];

%for shifting and averaging:
FRAMES = -Max_Num_frames:Max_Num_frames;
FRAMES = FRAMES';
j = 1;
%%
for i=1:Num_of_graph_files
    %prepare for shifting:
    RED = zeros(length(FRAMES),1);
    GREEN = zeros(length(FRAMES),1);
    %get graph file name:
    graph_file_name = List_of_graph_files(i).name;
    %get track id:
    pos_of_dot = strfind(convertCharsToStrings(graph_file_name),'.');
    track_id = str2double(graph_file_name(Len_graph_prefix+1:pos_of_dot-1));
    
    %get red and green tracks under a specific track id:
    track_id_row = find(tracks.id==track_id);%find the row in tracks containing track_id
    frames = tracks.frames{track_id_row};
    red = tracks.ch1int{track_id_row};
    green = tracks.ch2int{track_id_row};
    
    %find Red drop point:
    max_g_red = findDropPoint(green,red,pass_band);
    
    id = [id;track_id];

    if ~isempty(max_g_red)%if there is a Drop point
        id1 = [id1;track_id];
        drop_frame = [drop_frame;frames(max_g_red)];
        %shift the CrossOver point to 0:
        RED(ceil(length(FRAMES)/2)-max_g_red+1:ceil(length(FRAMES)/2)+length(red)-max_g_red,1) = red;
        GREEN(ceil(length(FRAMES)/2)-max_g_red+1:ceil(length(FRAMES)/2)+length(red)-max_g_red,1) = green;
        TRACKS{j,1} = RED;
        TRACKS{j,2} = GREEN;
        j = j + 1;
    else
        drop_frame = [drop_frame;-1000];
    end
end

%%
T = table(id,drop_frame);
T1 = sortrows(T);
%SaveNameE = [SaveName,'.csv'];
%writetable(T1,SaveNameE)
%disp(['Result file written to .csv in: ', pwd]);
%shifted tracks
[trackid,CH1INT,CH2INT]=sortcell(id1,TRACKS,FRAMES);
%normalization:
for i=1:length(CH1INT)
    for j=1:length(FRAMES)
        r(i,j) = CH1INT{i,1}(j);
        g(i,j) = CH2INT{i,1}(j);
    end
%     if max(r(i,:))>=max(g(i,:))
%         M = max(r(i,:));
%     else
%         M = max(g(i,:));
%     end
    pos = ceil(length(FRAMES)/2)+1;
    M = g(i,pos);
    R(i,:) = r(i,:)/M;%max(r(i,:));
    G(i,:) = g(i,:)/M;%max(g(i,:));
end
%find mean and median:
%mean_R = mean(R,1);
MEAN_R = MEAN(R);
MEDIAN_R = MEDIAN(R);
%mean_G = mean(G,1);
MEAN_G = MEAN(G);
MEDIAN_G = MEDIAN(G);
FRAMES = FRAMES/6;
c3 = findXticks(FRAMES,x_ticks);
x = (-c3*x_ticks:x_ticks:c3*x_ticks);

%plot red without removing any data point
FigH = figure('Position', get(0, 'Screensize'));
hold on
for i=1:size(R,1)
    xplot = FRAMES;
    yplot = R(i,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win2,'moving'),'color', [240/255,128/255,128/255])
end
%plot(FRAMES,MEAN_R,'k','LineWidth',2)
MEDIAN_R(MEDIAN_R==0) = NaN;
plot(FRAMES,smooth(MEDIAN_R,win2,'moving'),'color',[220/255,20/255,60/255],'LineWidth',3)
% plot(FRAMES,MEDIAN_R,'k','LineWidth',5)
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
% saveas(FigH,'RedChannel.fig')
hold off;
% close all;


R2 = R(:,ceil(length(FRAMES)/2):end);
G2 = G(:,ceil(length(FRAMES)/2):end);

clear temp;
temp = R2;
temp(temp == 0) = nan;
R3 = zeros(size(R2));
for i = 1:size(temp,1)
    R3(i,:) = smooth(temp(i,:),win2,'moving');
%     temp2 = (temp(~isnan(temp(i,:))));
%     temp2(isnan(temp2))=0;
%     R3(i,1:size(temp2,2)) = lowpass(temp2,0.1);
end
R3(R3==0) = NaN;
clear temp temp2;
temp = G2;
temp(temp == 0) = nan;
G3 = zeros(size(G2));
for i = 1:size(temp,1)
    G3(i,:) = smooth(temp(i,:),win2,'moving');
%     temp2 = (temp(~isnan(temp(i,:))));
%     temp2(isnan(temp2))=0;
%     G3(i,1:size(temp2,2)) = lowpass(temp2,0.1);
end
G3(G3==0) = NaN;
clear temp temp2;

MEAN_R2 = MEAN(R2);
MEDIAN_R2 = MEDIAN(R2);
%mean_G = mean(G,1);
MEAN_G2 = MEAN(G2);
MEDIAN_G2 = MEDIAN(G2);


RR = rmoutliers(R3,'percentiles',[low high]);
GG = rmoutliers(G3,'percentiles',[low high]);

RR(isnan(RR)) = 0;
GG(isnan(GG)) = 0;

MEAN_R3 = MEAN(RR);
MEDIAN_R3 = MEDIAN(RR);
%mean_G = mean(G,1);
MEAN_G3 = MEAN(GG);
MEDIAN_G3 = MEDIAN(GG);

FRAMES2 = FRAMES(ceil(length(FRAMES)/2):end);

%plot red where we remove points before zero and smooth the remaning curves
%before plotting

FigH = figure('Position', get(0, 'Screensize'));
hold on
for i=1:size(R2,1)
    xplot = FRAMES2;
    yplot = R2(i,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win2,'moving'),'color', [240/255,128/255,128/255])
end
%plot(FRAMES,MEAN_R,'k','LineWidth',2)
% MEDIAN_R2(MEDIAN_R2==0) = NaN;
% plot(FRAMES2,smooth(MEDIAN_R2,win2,'moving'),'color',[220/255,20/255,60/255],'LineWidth',3)
MEDIAN_R3(MEDIAN_R3==0) = NaN;
plot(FRAMES2,smooth(MEDIAN_R3,win2,'moving'),'color',[220/255,20/255,60/255],'LineWidth',3)
% plot(FRAMES,MEDIAN_R,'k','LineWidth',5)
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
% saveas(FigH,'RedChannel.fig')
hold off;

%plot red where we remove points before zero and smooth the remaning curves
%before plotting

% R3 = rmoutliers(R2,'percentiles',[5 95]);
% G3 = rmoutliers(G2,'percentiles',[5 95]);

% plot red removing the points before zero, smoothing and shaded error bars
FigH = figure('Position', get(0, 'Screensize'));
hold on;
shadedErrorBar(FRAMES2,rmoutliers(R3,'percentiles',[low high]),{@median,@std},'lineprops',{'-r', 'MarkerFaceColor', 'r'},'patchSaturation',patchSatVal);
% shadedErrorBar(FRAMES2,R2,{@median,@std},'lineprops',{'-r', 'MarkerFaceColor', 'r'},'patchSaturation',patchSatVal);

xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
% saveas(FigH,'RedChannel.fig')
hold off;


%plot green without removing any data point
FigH = figure('Position', get(0, 'Screensize'));
hold on
for i=1:size(G,1)
    xplot = FRAMES;
    yplot = G(i,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win2,'moving'),'color', [144/255,238/255,144/255])
end
%plot(FRAMES,MEAN_R,'k','LineWidth',2)
MEDIAN_G(MEDIAN_G==0) = NaN;
plot(FRAMES,smooth(MEDIAN_G,win2,'moving'),'color',[0/255, 128/255, 0/255],'LineWidth',3)
% plot(FRAMES,MEDIAN_G,'k','LineWidth',5)
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
% saveas(FigH,'GreenChannel.fig')
hold off;
% close all;


%plot green where we remove points before zero and smooth the remaning curves
%before plotting

FigH = figure('Position', get(0, 'Screensize'));
hold on
for i=1:size(G2,1)
    xplot = FRAMES2;
    yplot = G2(i,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win2,'moving'),'color', [144/255,238/255,144/255])
end
%plot(FRAMES,MEAN_R,'k','LineWidth',2)
% MEDIAN_G2(MEDIAN_G2==0) = NaN;
% plot(FRAMES2,smooth(MEDIAN_G2,win2,'moving'),'color',[0/255, 128/255, 0/255],'LineWidth',3)
MEDIAN_G3(MEDIAN_G3==0) = NaN;
plot(FRAMES2,smooth(MEDIAN_G3,win2,'moving'),'color',[0/255, 128/255, 0/255],'LineWidth',3)
% plot(FRAMES,MEDIAN_R,'k','LineWidth',5)
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
% saveas(FigH,'RedChannel.fig')
hold off;

% plot green removing the points before zero, smoothing and shaded error bars
FigH = figure('Position', get(0, 'Screensize'));
hold on;
shadedErrorBar(FRAMES2,rmoutliers(G3,'percentiles',[low high]),{@median,@std},'lineprops',{'-g', 'MarkerFaceColor', 'g'},'patchSaturation',patchSatVal);
% shadedErrorBar(FRAMES2,G2,{@median,@std},'lineprops',{'-g', 'MarkerFaceColor', 'g'},'patchSaturation',patchSatVal);

xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
% saveas(FigH,'RedChannel.fig')
hold off;

% Plot Green and Red removing points before zero, including their shaded
% error bars

%plot(FRAMES,mean_R,'r',FRAMES,mean_G,'g')
%figure,plot(FRAMES,MEAN_R,'r',FRAMES,MEAN_G,'g')
FigH = figure('Position', get(0, 'Screensize'));
% plot(FRAMES,MEDIAN_R,'r',FRAMES,MEDIAN_G,'g')
hold on;
shadedErrorBar(FRAMES2,rmoutliers(G3,'percentiles',[low high]),{@median,@std},'lineprops',{'-g', 'MarkerFaceColor', 'g'},'patchSaturation',patchSatVal);
shadedErrorBar(FRAMES2,rmoutliers(R3,'percentiles',[low high]),{@median,@std},'lineprops',{'-r', 'MarkerFaceColor', 'r'},'patchSaturation',patchSatVal);
% shadedErrorBar(FRAMES2,G2,{@median,@std},'lineprops',{'-g', 'MarkerFaceColor', 'g'},'patchSaturation',patchSatVal);
% shadedErrorBar(FRAMES2,R2,{@median,@std},'lineprops',{'-r', 'MarkerFaceColor', 'r'},'patchSaturation',patchSatVal);
% shadedErrorBar(FRAMES2,lowpass(G3,0.05),{@median,@std},'lineprops',{'-g', 'MarkerFaceColor', 'g'},'patchSaturation',patchSatVal);
% shadedErrorBar(FRAMES2,lowpass(R3,0.05),{@median,@std},'lineprops',{'-r', 'MarkerFaceColor', 'r'},'patchSaturation',patchSatVal);

xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
%title('Graph of Mean Value')
title('Graph of Median Value')
% saveas(FigH,SaveGraphName)
hold off;
save ([csv_path,'All.mat']);
%%
rmpath(graph.dir);
toc

%%
function max_g_red = findDropPoint(green,red,pass_band)
    %lowpass filtering of tracks:
    f_green = lowpass(green,pass_band);
    f_red = lowpass(red,pass_band);
    
    %find gradients of lowpass filtered tracks:
    g_green = gradient(f_green);
    g_red = gradient(f_red);
    gg_red = gradient(g_red);
     
    %find the maximum of 2nd derivative of the lowpass filtered red track:
    if((max(gg_red)>0.1) && (g_red(gg_red==max(gg_red)) < 0))
        max_g_red = find(gg_red==max(gg_red));
    else
        max_g_red = find(gg_red~=0,1,'first');
    end
    
%     %%find the right frame where crossover occurs
%     TF = 1;
%     while TF
%         %find frames with minus gradients:
%         frame_of_minus_gradient = find(g_green<0);%
%         if isempty(frame_of_minus_gradient)
%             min_g_green = [];
%             break;
%         end
%         %find the consecutive sections in frames with minus gradients:
%         [lenconsecuframe,endconsecuframe] = findconsecvec(frame_of_minus_gradient');
%         %find pos with minimum minus gradients in frame_of_minus_gradient:
%         pos_min_g_green = find(frame_of_minus_gradient==find(g_green==g_green(min_g_green)));
%         %find which consecutive section min_g_green belongs to: 
%         a1 = find(endconsecuframe>=find(frame_of_minus_gradient==find(g_green==g_green(min_g_green))));
%         startframe = endconsecuframe(a1(1)) - lenconsecuframe(a1(1)) + 1;
%         %find the frame interval within which crossover may occur:
%         frame_interval_of_interest = frame_of_minus_gradient(startframe:endconsecuframe(a1(1)),1);
%         
%         %if Red decreases and Green increases in the frame_interval_of_interest
%         if red(frame_interval_of_interest(end))<green(frame_interval_of_interest(end))||red(frame_interval_of_interest(1))>green(frame_interval_of_interest(1))
%             g_green(frame_interval_of_interest) = 0;
%             if ~isempty(find(g_green<0,1))
%                 min_g_green = find(g_green==min(g_green));
%             end
%         else
%             TF = 0;
%         end
%     end
end

function [c,d] = findconsecvec(q)
a=diff(q);
b=find([a inf]>1);%
c=diff([0 b]); %length of the sequences
d=cumsum(c); %endpoints of the sequences
end

function [Max_Num_frames,List_of_graph_files] = findMaxNumFrames(graph,tracks)
Max_Num_frames = 0;
Len_graph_prefix = length(graph.prefix);
List_of_graph_files = dir([graph.dir,'*.',graph.suffix]);%
Num_of_graph_files = length(List_of_graph_files);%number of graph files
for i=1:Num_of_graph_files
    %get graph file name:
    graph_file_name = List_of_graph_files(i).name;
    %get track id:
    pos_of_dot = strfind(convertCharsToStrings(graph_file_name),'.');
    track_id = str2double(graph_file_name(Len_graph_prefix+1:pos_of_dot-1));
    %find the row in tracks containing track_id:
    track_id_row = find(tracks.id==track_id);
    %pick frame number series with a track defined by track_id:
    frames = tracks.frames{track_id_row};
    %find the Maximum number of framees among all tracks define by all track_id:
    if Max_Num_frames<=length(frames)%find max num of frames among all tracks
        Max_Num_frames = length(frames);
    end
end
end

function [id1,CH1INT,CH2INT] = sortcell(id,TRACKS,FRAMES)
%Sort a MxN cell array 'TRACKS' according to a double array 'id'
i = 1;
while ~isempty(id)
    min_id = min(id);
    min_id_idx = find(id==min_id);
    id1(i) = id(min_id_idx);
    id(min_id_idx) = [];
    CH1INT{i,1} = TRACKS{min_id_idx,1};
    CH2INT{i,1} = TRACKS{min_id_idx,2};
    TRACKS(min_id_idx,:) = [];
    i = i + 1;
end
end

function MEAN_R = MEAN(R)
sR = size(R);%86,471
for j=1:sR(2)
    a = R(:,j);
    an = nonzeros(a);
    if ~isempty(an)
        MEAN_R(j) = mean(an);
    else
        MEAN_R(j) = 0;
    end
end
end

function MEDIAN_R = MEDIAN(R)
sR = size(R);%86,471
for j=1:sR(2)
    a = R(:,j);
    an = nonzeros(a);
    if ~isempty(an)
        MEDIAN_R(j) = median(an);
    else
        MEDIAN_R(j) = 0;
    end
end
end

function c3 = findXticks(FRAMES,x_ticks)
c = (length(FRAMES) - 1)/2;
c1 = 0;
c2 = -1;
while c2<0
c1 = c1 + 1;
c2 = c1*x_ticks - c;
end
c3 = c1 - 1;
end