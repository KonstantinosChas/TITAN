% This is the main routine to handle data from wave gauges, wave paddles, 
% microsonic distance sensor and hydrophone.
% For each instrument, data are being loaded by sub-routines (.m) 
% which are then being nested in analysis/plotting - functions
% These parent functions are being called here. These function load data, 
% do some preliminary analysis and plot the basics
% Add to the current script for comparison between instruments;
% For further analysis add to the existing functions or create new ones 
% and then call them from here.


% To use the current script properly you should place the folder-data files for each instruments 
% in the appropriate folders designated in the paths below e.g. [mainpath,instrument_path].

%% Section 1: paths 
clear; clc;
scrp_pth = 'C:\Users\kchasapi\Desktop\konstantinos\scripts\matlab';
addpath(scrp_pth);
main_pth = 'C:\Users\kchasapi\Desktop\konstantinos\data';
pd_pth = '\surf_elev_pd\';
wg_pth = '\surf_elev\';
micro_pth = '\micro\';
hydro_pth = '\hydro\';

dir_plot_save_pd = ['C:\Users\kchasapi\Desktop\konstantinos\plots',pd_pth];
dir_plot_save_wg = ['C:\Users\kchasapi\Desktop\konstantinos\plots',wg_pth];
dir_plot_save_mi = ['C:\Users\kchasapi\Desktop\konstantinos\plots',micro_pth];
dir_plot_save_hy = ['C:\Users\kchasapi\Desktop\konstantinos\plots',hydro_pth];


%% Section 2: calling all instruments
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % wave paddles 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% listdir = ('*.txt');
% numwg = length(listdir);
% header = zeros(1,numwg) + 5;
% wgnum = zeros(1,numwg) + 16; 
% sf = zeros(1,numwg) + 128; 
% junktime = zeros(1,numwg) + 4;
% tottime = zeros(1,numwg) + 64;
% for i = numwg
%     pd_full= fullfile(main_pth,pd_pth);
% end
% [pd_data] = read_analyse_pddle_dwbasin(pd_full,dir_plot_save_pd,header,tottime,sf);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % surface elevation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = numwg
%     wg_full= fullfile(main_pth,wg_pth);
% end
% [wg_data] = read_analyse_wg_dwbasin(wg_full,dir_plot_save_wg,header,wgnum,tottime, sf,junktime);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % microsonic 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% micro_fld = {'mi_20_07_16_13_45','mi_20_07_16_46_10','mi_20_07_17_18_35','mi_20_07_17_54_46','mi_20_07_18_27_11','mi_20_07_18_58_55'};
% nummicro = length(micro_fld);
% micro_full = cell(1,nummicro);
% for i = 1:nummicro
%     micro_full{i} = fullfile(main_pth,micro_pth,micro_fld{i});
% end
% [micro_data] = read_analyse_micro(micro_full,dir_plot_save_mi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hydrophone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hydro_fld = {'hy_20_07_16_14_05','hy_20_07_16_46_30','hy_20_07_17_18_14'};
numhydro = length(hydro_fld);
hydro_full = cell(1,numhydro);
for i = 1:numhydro
    hydro_full{i} = fullfile(main_pth,hydro_pth,hydro_fld{i});
end
[hydro_data] = read_analyse_hydro(hydro_full,dir_plot_save_hy);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 3: comparisons between instruments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







