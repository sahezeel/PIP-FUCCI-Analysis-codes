function FindCrossOver(csv_path)
%subfolder containing graph files, this subfolder should be in the same folder as this matlab file
graph.subdir = 'good G1';
%prefix and suffix of graph file names:
graph.prefix = 'Track ID';
graph.suffix = 'png';

%save names:
SaveGraphName = 'MedianGraph.fig';
%SaveName = 'CrossOverFrames1';

%pass band for lowpass filter:
pass_band = 0.1;%choosing from (0,1)

%xticks
x_ticks = 2;

%window size for moving average filter to smooth the graphs
win = 7; %individual graph smoothning
win2 = 7; % median smoothening

%The color of the confidence intervals is decided by setting this value
%between 0 and 1. 0 is blank with no color and 1 is completely rich
patchSatVal = 0.1;
M = 6;

%%
%Preparation:
graph.dir = [csv_path,graph.subdir,'\'];
addpath(graph.dir);
Len_graph_prefix = length(graph.prefix);%length of prefix of graph file names
tracks_file = [csv_path,'tracks.mat'];
load(tracks_file,'tracks')% tracks;%tracks

%%
%list graph files:
%find max num of frames in all tracks defined by graph files:
[Max_Num_frames,List_of_graph_files] = findMaxNumFrames(graph,tracks);
Num_of_graph_files = length(List_of_graph_files);%number of graph files

%for crossover:
id = [];
id1 = [];
cross_over_frame = [];

%for shifting and averaging:
FRAMES = -Max_Num_frames:Max_Num_frames;
FRAMES = FRAMES';
j1 = 1;
%%
for i1=1:Num_of_graph_files
    %prepare for shifting:
    RED = zeros(length(FRAMES),1);
    GREEN = zeros(length(FRAMES),1);
    %get graph file name:
    graph_file_name = List_of_graph_files(i1).name;
    %get track id:
    pos_of_dot = strfind(convertCharsToStrings(graph_file_name),'.');
    track_id = str2double(graph_file_name(Len_graph_prefix+1:pos_of_dot-1));
    
    %get red and green tracks under a specific track id:
    track_id_row = find(tracks.id==track_id);%find the row in tracks containing track_id
    frames = tracks.frames{track_id_row};
    red = tracks.ch1int{track_id_row};
    green = tracks.ch2int{track_id_row};
    
    %find CrossOver point:
    min_g_green = findCrossOver(green,red,pass_band);
    
    id = [id;track_id];

    if ~isempty(min_g_green)%if there is a CrossOver point
        id1 = [id1;track_id];
        cross_over_frame = [cross_over_frame;frames(min_g_green)];
        %shift the CrossOver point to 0:
        RED(ceil(length(FRAMES)/2)-min_g_green+1:ceil(length(FRAMES)/2)+length(red)-min_g_green,1) = red;
        GREEN(ceil(length(FRAMES)/2)-min_g_green+1:ceil(length(FRAMES)/2)+length(red)-min_g_green,1) = green;
        TRACKS{j1,1} = RED;
        TRACKS{j1,2} = GREEN;
        j1 = j1 + 1;
    else
        cross_over_frame = [cross_over_frame;-1000];
    end
end

%%
T = table(id,cross_over_frame);
T1 = sortrows(T);
%SaveNameE = [SaveName,'.csv'];
%writetable(T1,SaveNameE)
%disp(['Result file written to .csv in: ', pwd]);
%shifted tracks
[trackid,CH1INT,CH2INT]=sortcell(id1,TRACKS,FRAMES);
%normalization:

%---------------- You can comment this part -----------------------------%
for i=1:length(CH1INT)
    for j=1:length(FRAMES)
        r(i,j) = CH1INT{i,1}(j);
        g(i,j) = CH2INT{i,1}(j);
    end
    if max(r(i,:))>=max(g(i,:))
        M = max(r(i,:));
    else
        M = max(g(i,:));
    end
    R(i,:) = r(i,:)/M;%max(r(i,:));
    G(i,:) = g(i,:)/M;%max(g(i,:));
end
%-------------------------------------------------------------------------%
%--------------- And Uncomment this one ----------------------------------%
% for i=1:length(CH1INT)
%     for j=1:length(FRAMES)
%         r(i,j) = CH1INT{i,1}(j);
%         g(i,j) = CH2INT{i,1}(j);
%     end
%     if max(r(i,ceil(length(FRAMES)/2)-6:ceil(length(FRAMES)/2)))>=max(g(i,ceil(length(FRAMES)/2)-6:ceil(length(FRAMES)/2)))
%         M = max(r(i,ceil(length(FRAMES)/2)-6:ceil(length(FRAMES)/2)));
%     else
%         M = max(g(i,ceil(length(FRAMES)/2)-6:ceil(length(FRAMES)/2)));
%     end
%     R(i,:) = r(i,:)/M;%max(r(i,:));
%     G(i,:) = g(i,:)/M;%max(g(i,:));
% end
%-------------------------------------------------------------------------%
clear temp;
temp = G;
temp(temp == 0) = nan;
for i = 1:size(temp,1)
    G2(i,:) = smooth(temp(i,:),win,'moving');
end
clear temp;
temp = R;
temp(temp == 0) = nan;
for i = 1:size(temp,1)
    R2(i,:) = smooth(temp(i,:),win,'moving');
end
clear temp;

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

%plot red with outliers 

FigH = figure('Position', get(0, 'Screensize'));
hold on
for i2=1:size(R,1)
    xplot = FRAMES;
    yplot = R(i2,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win,'moving'),'color', [220/255,20/255,60/255],'LineWidth',1) %tomato individual
end
% plot(FRAMES,smooth(MEAN_R,win2,'moving'),'color',[220/255,20/255,60/255],'LineWidth',2)
plot(FRAMES,smooth(MEDIAN_R,win2,'moving'),'color',[220/255,20/255,60/255],'LineWidth',6) %red original
xticks(x)
xlabel('Time (Hours)','FontSize',12,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',12,'FontWeight','bold')
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
savefile = [csv_path,'RedChannel_w_outliers.fig'];
%saveas(FigH,savefile)
savefig(FigH,savefile)
hold off

% Plot Red channel removing the outliers

R3 = rmoutliers(R2,'percentiles',[1 99]);
FigH = figure('Position', get(0, 'Screensize'));
hold on
for i=1:size(R3,1)
    xplot = FRAMES;
    yplot = R3(i,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win,'moving'),'color', [220/255,20/255,60/255])
end
% plot(FRAMES,smooth(MEAN_R,win2,'moving'),'color',[220/255,20/255,60/255],'LineWidth',2)
plot(FRAMES,smooth(MEDIAN_R,win2,'moving'),'color',[220/255,20/255,60/255],'LineWidth',6)
xticks(x)
xlabel('Time (Hours)','FontSize',12,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',12,'FontWeight','bold')
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
savefile = [csv_path,'RedChannel_wo_outliers.fig'];
%saveas(FigH,savefile)
savefig(FigH,savefile)
hold off

% Plot red without outliers and shaded error bars

FigH = figure('Position', get(0, 'Screensize'));
hold on;
shadedErrorBar(FRAMES,rmoutliers(R2,'percentiles',[1 99]),{@median,@std},'lineprops',{'-','color', [220/255,20/255,60/255], 'MarkerFaceColor', 'r', 'LineWidth', M},'patchSaturation',patchSatVal);
xticks(x)
xlabel('Time (Hours)','FontSize',12,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',12,'FontWeight','bold')
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
savefile = [csv_path,'RedChannel_shaded.fig'];
%saveas(FigH,savefile)
savefig(FigH,savefile)
hold off


%plot green with outliers 

FigH = figure('Position', get(0, 'Screensize'));
hold on
for i2=1:size(G,1)
    xplot = FRAMES;
    yplot = G(i2,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win,'moving'),'color', [34/255,139/255,34/255] ,'LineWidth',1)
end
% plot(FRAMES,smooth(MEAN_G,win2,'moving'),'color',[0/255, 128/255, 0/255],'LineWidth',2)
plot(FRAMES,smooth(MEDIAN_G,win2,'moving'),'color',[0/255, 100/255, 0/255],'LineWidth',6)
xticks(x)
xlabel('Time (Hours)','FontSize',12,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',12,'FontWeight','bold')
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
savefile = [csv_path,'GreenChannel_w_outliers.fig'];
%saveas(FigH,savefile)
savefig(FigH,savefile)
hold off

% Plot Green channel removing the outliers

G3 = rmoutliers(G2,'percentiles',[1 99]);
FigH = figure('Position', get(0, 'Screensize'));
hold on
for i=1:size(G3,1)
    xplot = FRAMES;
    yplot = G3(i,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win,'moving'),'color', [144/255,238/255,144/255]) %individual plots change colors here
end
% plot(FRAMES,smooth(MEAN_G,win2,'moving'),'color',[0/255, 128/255, 0/255],'LineWidth',2)
plot(FRAMES,smooth(MEDIAN_G,win2,'moving'),'color',[0/255, 128/255, 0/255],'LineWidth',6)
xticks(x)
xlabel('Time (Hours)','FontSize',12,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',12,'FontWeight','bold')
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
savefile = [csv_path,'GreenChannel_wo_outliers.fig'];
%saveas(FigH,savefile)
savefig(FigH,savefile)
hold off

% Plot green without outliers and shaded error bars

FigH = figure('Position', get(0, 'Screensize'));
hold on;
shadedErrorBar(FRAMES,rmoutliers(G2,'percentiles',[1 99]),{@median,@std},'lineprops',{'-', 'color' [0/255,128/255,0/255] 'MarkerFaceColor', 'g', 'LineWidth', M},'patchSaturation',patchSatVal);
xticks(x)
xlabel('Time (Hours)','FontSize',12,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',12,'FontWeight','bold')
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
savefile = [csv_path,'GreenChannel_shaded.fig'];
%saveas(FigH,savefile)
savefig(FigH,savefile)
hold off


%plot(FRAMES,mean_R,'r',FRAMES,mean_G,'g')
%figure,plot(FRAMES,MEAN_R,'r',FRAMES,MEAN_G,'g')
FigH = figure('Position', get(0, 'Screensize'));
hold on;
shadedErrorBar(FRAMES,rmoutliers(G2,'percentiles',[1 99]),{@median,@std},'lineprops',{'-', 'color', [0/255,128/255,0/255], 'MarkerFaceColor', 'g', 'LineWidth', M},'patchSaturation',patchSatVal);

shadedErrorBar(FRAMES,rmoutliers(R2,'percentiles',[1 99]),{@median,@std},'lineprops',{'-', 'color', [220/255,20/255,60/255], 'MarkerFaceColor', 'r', 'LineWidth', M},'patchSaturation',patchSatVal);
xticks(x)
xlabel('Time (Hours)','FontSize',12,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',12,'FontWeight','bold')
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
set(gca,'XTickLabelMode','auto')
%title('Graph of Mean Value')
title('Graph of Median Value','FontSize',50,'FontWeight','bold')
savefile = [csv_path,SaveGraphName];
%saveas(FigH,savefile)
savefig(FigH,savefile)
% save([csv_path,'MEANS.mat'],'FRAMES','MEAN_R','MEAN_G','x')
save([csv_path,'MEDIANS.mat'],'FRAMES','MEDIAN_R','MEDIAN_G','x', 'G', 'G3', 'R', 'R3');

%%
rmpath(graph.dir);
close all
save ([csv_path,'Entire workspace.mat']);

%%
function min_g_green = findCrossOver(green,red,pass_band)
    %lowpass filtering of tracks:
    f_green = lowpass(green,pass_band);
    f_red = lowpass(red,pass_band);
    
    %find gradients of lowpass filtered tracks:
    g_green = gradient(f_green);
    g_red = gradient(f_red);
     
    %find the minimum of gradients of lowpass filtered tracks:
    min_g_green = find(g_green==min(g_green));
    
    %%find the right frame where crossover occurs
    TF = 1;
    while TF
        %find frames with minus gradients:
        frame_of_minus_gradient = find(g_green<0);%
        if isempty(frame_of_minus_gradient)
            min_g_green = [];
            break;
        end
        %find the consecutive sections in frames with minus gradients:
        [lenconsecuframe,endconsecuframe] = findconsecvec(frame_of_minus_gradient');
        %find pos with minimum minus gradients in frame_of_minus_gradient:
        pos_min_g_green = find(frame_of_minus_gradient==find(g_green==g_green(min_g_green)));
        %find which consecutive section min_g_green belongs to: 
        a1 = find(endconsecuframe>=find(frame_of_minus_gradient==find(g_green==g_green(min_g_green))));
        startframe = endconsecuframe(a1(1)) - lenconsecuframe(a1(1)) + 1;
        %find the frame interval within which crossover may occur:
        frame_interval_of_interest = frame_of_minus_gradient(startframe:endconsecuframe(a1(1)),1);
        
        %if Red decreases and Green increases in the frame_interval_of_interest
        if red(frame_interval_of_interest(end))<green(frame_interval_of_interest(end))||red(frame_interval_of_interest(1))>green(frame_interval_of_interest(1))
            g_green(frame_interval_of_interest) = 0;
            if ~isempty(find(g_green<0,1))
                min_g_green = find(g_green==min(g_green));
            end
        else
            TF = 0;
        end
    end
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
for i4=1:Num_of_graph_files
    %get graph file name:
    graph_file_name = List_of_graph_files(i4).name;
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
i5 = 1;
while ~isempty(id)
    min_id = min(id);
    min_id_idx = find(id==min_id);
    id1(i5) = id(min_id_idx);
    id(min_id_idx) = [];
    CH1INT{i5,1} = TRACKS{min_id_idx,1};
    CH2INT{i5,1} = TRACKS{min_id_idx,2};
    TRACKS(min_id_idx,:) = [];
    i5 = i5 + 1;
end
end

function MEAN_R = MEAN(R)
sR = size(R);%86,471
for j2=1:sR(2)
    a = R(:,j2);
    an = nonzeros(a);
    if ~isempty(an)
        MEAN_R(j2) = mean(an);
    else
        MEAN_R(j2) = 0;
    end
end
end

function MEDIAN_R = MEDIAN(R)
sR = size(R);%86,471
for j3=1:sR(2)
    a = R(:,j3);
    an = nonzeros(a);
    if ~isempty(an)
        MEDIAN_R(j3) = median(an);
    else
        MEDIAN_R(j3) = 0;
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

%%
end