clear all
close all
clc

%% Input data for the simulink simulation
table = readtable("./../data/reference.tsv", "FileType","text",'Delimiter', '\t');

simin = table2array([table(:,1) table(:,2)]);

% Now run the simulink model!