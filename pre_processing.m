clear
clc
close all

%%
%==========================User Inputs==============================
%path to the .csv file:
csv_path = '\\corefs2.med.umich.edu\Shared3\ALR_Lab\SA_FADD_Data\2020-2021\1403 ZEN Confocal DATA ALL\08272021 A2 FADDKD Palblo MTX Tram test\Analysis plot raw green and red\16\';
csv_file = 'All Spots statistics.csv';

%%
path = [csv_path,csv_file];
%create graph folder:
figpath = [csv_path,'graph_200\'];
if exist(figpath,'dir')
    rmdir(figpath,'s')
end
mkdir(figpath)
good_G1_graph = [csv_path,'good G1\'];
if ~exist(good_G1_graph,'dir')
    mkdir(good_G1_graph)
end
%
tracks_file = [csv_path,'tracks.mat'];
if isfile(tracks_file)
    delete(tracks_file)
end

%save tracks:
ParseSpotStats2;
save(tracks_file,'tracks')
%save png graphs:
FilterIntensities1;