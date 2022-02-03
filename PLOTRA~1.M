
function PlotRawGreenRed(csv_path)
% close all 
% clear, clc 
% tic
%
%%===================   ===== User inputs ============================
%path to the .csv file:
% csv_path = 'C:\Users\sundarer\Desktop\GalbanLab\Sahezeel_Al_MIL\New Test Data for Ram\Analysis using Zhouping codes\Normal cell cycle\H1975 CT 16\';

%subfolder containing graph files, this subfolder should be in the same folder as this matlab file
graph.subdir = 'graph';
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
win = 50;
win2 = 10;
win3 = 100;
win4 = 250;

%The color of the confidence intervals is decided by setting this value
%between 0 and 1. 0 is blank with no color and 1 is completely rich
patchSatVal = 0.1;

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
drop_frame = [];

%for shifting and averaging:
% FRAMES = -Max_Num_frames:Max_Num_frames;
FRAMES = 0:Max_Num_frames;
FRAMES = FRAMES';
j = 1;
%%
t = [];
for i=1:Num_of_graph_files
    %prepare for shifting:
    RED = zeros((length(FRAMES)+300),1);
    GREEN = zeros((length(FRAMES)+300),1);
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
%     RED = tracks.ch1int{track_id_row};
%     GREEN = tracks.ch2int{track_id_row};
    t = [t; length(frames)];
    
    %find Red drop point:
%     max_g_red = findDropPoint(green,red,pass_band);
    
%     id = [id;track_id];

%     if ~isempty(max_g_red)%if there is a Drop point
    id1 = [id1;track_id];
%     drop_frame = [drop_frame;frames(max_g_red)];
        %shift the CrossOver point to 0:
%     RED(ceil(length(FRAMES)/2)-max_g_red+1:ceil(length(FRAMES)/2)+length(red)-max_g_red,1) = red;
%     GREEN(ceil(length(FRAMES)/2)-max_g_red+1:ceil(length(FRAMES)/2)+length(red)-max_g_red,1) = green;
    RED(frames(1)+1:frames(1)+length(frames)) = red;
    GREEN(frames(1)+1:frames(1)+length(frames)) = green;
    TRACKS{j,1} = RED;
    TRACKS{j,2} = GREEN;
    j = j + 1;
%     else
%         drop_frame = [drop_frame;-1000];
%     end
end

%%
% T = table(id,drop_frame);
% T1 = sortrows(T);
%SaveNameE = [SaveName,'.csv'];
%writetable(T1,SaveNameE)
%disp(['Result file written to .csv in: ', pwd]);
%shifted tracks
[trackid,CH1INT,CH2INT]=sortcell(id1,TRACKS,FRAMES);
%normalization:
r = zeros(length(CH1INT),length(FRAMES));
g = zeros(length(CH1INT),length(FRAMES));
for i=1:length(CH1INT)
    for j=1:length(CH1INT{i,1})
        r(i,j) = CH1INT{i,1}(j);
        g(i,j) = CH2INT{i,1}(j);
    end
%     if max(r(i,:))>=max(g(i,:))
%         M = max(r(i,:));
%     else
%         M = max(g(i,:));
%     end
%     pos = ceil(length(FRAMES)/2)+1;
%     M = g(i,pos);
    M = 1;
    R(i,:) = r(i,:)/M;%max(r(i,:));
    G(i,:) = g(i,:)/M;%max(g(i,:));
end
%find mean and median:
%mean_R = mean(R,1);
FRAMES = 0:(Max_Num_frames+300);
MEAN_R = MEAN(R);
MEDIAN_R = MEDIAN(R);
%mean_G = mean(G,1);
MEAN_G = MEAN(G);
MEDIAN_G = MEDIAN(G);
FRAMES = FRAMES/6;
c3 = findXticks(FRAMES,x_ticks);
x = (-c3*x_ticks:x_ticks:c3*x_ticks);

%---------------plot red without removing any data point------------------%
FigH = figure('Position', get(0, 'Screensize'));
hold on
for i=1:size(R,1)
    xplot = FRAMES;
    yplot = R(i,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win,'moving'),'color', [240/255,128/255,128/255])
end
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
savefile = [csv_path,'RedChannel.fig'];
saveas(FigH,savefile)
hold off;
% close all;
%-------------------------------------------------------------------------%

%---------plot red without removing any data point with Median------------%
FigH = figure('Position', get(0, 'Screensize'));
hold on
for i=1:size(R,1)
    xplot = FRAMES;
    yplot = R(i,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win,'moving'),'color', [240/255,128/255,128/255])
end
%plot(FRAMES,MEAN_R,'k','LineWidth',2)
MEDIAN_R(MEDIAN_R==0) = NaN;
plot(FRAMES,smooth(MEDIAN_R,win2,'moving'),'color',[220/255,20/255,60/255],'LineWidth',3)
% plot(FRAMES,MEDIAN_R,'k','LineWidth',5)
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
savefile = [csv_path,'RedChannel_w_Median.fig'];
saveas(FigH,savefile)
hold off;
% close all;
%-------------------------------------------------------------------------%



%--------------plot red with median and shaded error bars-----------------%
FigH = figure('Position', get(0, 'Screensize'));
hold on;
shadedErrorBar(FRAMES,R,{@median,@std},'lineprops',{'-', 'MarkerFaceColor', 'r', 'color', [220/255,20/255,60/255]},'patchSaturation',patchSatVal);
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
savefile = [csv_path,'RedChannel_shaded.fig'];
saveas(FigH,savefile)
hold off;
%-------------------------------------------------------------------------%


%---------------plot green without removing any data point----------------%
FigH = figure('Position', get(0, 'Screensize'));
hold on
for i=1:size(G,1)
    xplot = FRAMES;
    yplot = G(i,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win,'moving'),'color', [144/255,238/255,144/255])
end
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
savefile = [csv_path,'GreenChannel.fig'];
saveas(FigH,savefile)
hold off;
% close all;
%-------------------------------------------------------------------------%

%--------plot green without removing any data point with Median-----------%
FigH = figure('Position', get(0, 'Screensize'));
hold on
for i=1:size(G,1)
    xplot = FRAMES;
    yplot = G(i,:);
    yplot(yplot==0) = NaN;
    plot(xplot,smooth(yplot,win,'moving'),'color', [144/255,238/255,144/255])
end
%plot(FRAMES,MEAN_R,'k','LineWidth',2)
MEDIAN_G(MEDIAN_G==0) = NaN;
plot(FRAMES,smooth(MEDIAN_G,win2,'moving'),'color',[0/255, 128/255, 0/255],'LineWidth',3)
% plot(FRAMES,MEDIAN_R,'k','LineWidth',5)
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
savefile = [csv_path,'GreenChannel_w_Median.fig'];
saveas(FigH,savefile)
hold off;
% close all;
%-------------------------------------------------------------------------%



%------------plot green with median and shaded error bars-----------------%
FigH = figure('Position', get(0, 'Screensize'));
hold on;
shadedErrorBar(FRAMES,G,{@median,@std},'lineprops',{'-', 'MarkerFaceColor', 'g', 'color', [0/255, 128/255, 0/255]},'patchSaturation',patchSatVal);
xticks(x)
xlabel('Time (Hours)','FontSize',20,'FontWeight','bold')
ylabel({'APC Activity';'(Normalized to max)'},'FontSize',20,'FontWeight','bold')
savefile = [csv_path,'GreenChannel_shaded.fig'];
saveas(FigH,savefile)
hold off;
%-------------------------------------------------------------------------%

save([csv_path,'MEDIANS_raw.mat'],'FRAMES','MEDIAN_R','MEDIAN_G','x', 'G', 'R');
%-------------------------------------------------------------------------%
%
rmpath(graph.dir);
% toc
close all;
save ([csv_path,'Entire workspace.mat']);

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
temp = [];
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
%     temp = [temp; length(frames)];
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

end