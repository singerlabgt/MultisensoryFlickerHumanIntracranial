%Produce all NatureComm2024 figures' panels.

%% Set dataset directory and dependent repos:

[~,repo]=define_flicker_root_dir;
for field_name=fieldnames(repo)'
    addpath(genpath(repo.(field_name{:})));
end

%note: can use function fetch_flicker_subjectIDs to get list of sessions
%and their details.

%%
clear all; close all; clc;
Figure_1_2_S2_S3A;

%%
clear all; close all; clc;
Figure_1_S2_S3A;

%%
clear all; close all; clc;
Figure_3_S3B;

%%
clear all; close all; clc;
plot_Figure_4A=0;
plot_Figure_4D_1=0;
Figure_4;

%%
clear all; close all; clc;
Figure_5;

%%
clear all; close all; clc;
plot_Figure_6A=0;
if plot_Figure_6A
    Figure_6A;
end

%%
clear all; close all; clc;
Figure_6B_S9;

%%
clear all; close all; clc;
Figure_S1;

%%
clear all; close all; clc;
Figure_S4;

%%
clear all; close all; clc;
Figure_S5;

%%
clear all; close all; clc;
Figure_S6;

%%
clear all; close all; clc;
Figure_S7;

%%
clear all; close all; clc;
Figure_S8;
