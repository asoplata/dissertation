%{
plotting_script: This script will recreate all of the subplots necessary
for Chapter 2 from scratch, so that they may be assembled using Inkscape.
Make sure to run "run_sims_and_calc_PAC_signal.m" first.

If you don't have data files of the correct filename in the current
directory, then this calls, for each missing data file, the individual
simulation-running function files to perform the actual simulations
and analysis; however, all of the plotting will happen here.

Dependencies: 
- specific DynaSim version installed (see README)
- "dynasim-extended-benita-model" model files for DynaSim installed (see README)
- A "recent" version of Matlab (both 2017a and 2018b were used originally)
- A "modern" desktop/laptop CPU (i.e. built in the 2010's)
- At LEAST 16GB of RAM (preferably more)
- At least 80GB of disk space
- Several hours of time to run the simulations (probably only 20 minutes 
each), run the analysis (long), and save the data (long). Depending on your
CPU and disk write speed, the entire data generation process may take 3 to 
6 hours or more. Of course, the time to generate each individual plot after loading
the proper data should take only a couple seconds.

%}

%% Part 0: Configuration

% Where to save all your simulation output subdirectories and data. This is
% set to the current directory by default. Note: Change the forward-slash
% to back-slash ('\') if on Windows.
output_dir = strcat(pwd, '/');

%% Part 1: Simulating and plotting wake/no-propofol case
 
% % The simulations can take a while to run (see above), and this script will
% % generate a LOT of figures, so it's recommended that you "step" through
% % the script manually if you want to really understand what's going on,
% % rather than running the whole script all at once. This is also very
% % useful for knowing what step you were on before Matlab inevitably
% % crashes, so you can continue where you left off. Relatedly: if you're
% % reading this, move away from Matlab to Python as soon as possible!
% % Comment out the line below if you DO want to run the whole script
% % non-stop:
% keyboard
% 
% Load the data if isn't present, and if the data file doesn't exist, run
% the simulation, its analyses, and save its data.
if ~exist('data_1_wake','var')
    if isfile('data_1_wake.mat')
        % Try to load the simulation data, if you've already run the sim
        % and it exists.
        load(strcat(output_dir,'data_1_wake.mat'))
    else
        % All code and configuration necessary to run and save both the
        % simulation and analyses is inside the following function file.
        % This will produce a folder containing around 17GB of data on
        % disk.
        data_1_wake = sim_1_wake_run(output_dir);
        % We'll clear our memory so we can combine and clean all the
        % saved-to-disk datasets
        clear data_1_wake

        % This will combine and clean the different data sets into a single
        % one, which we will save locally. After this finishes, you can
        % feel free to delete the 'sim_1_wake_run' folder containing
        % approximately 17GB of data sets. Note that this function
        % requires at LEAST 16GB of RAM available, so it's recommended to
        % close your other programs!
        data_1_wake = combine_data_files(output_dir, 'sim_1_wake');
        save(strcat(output_dir,'data_1_wake.mat'), 'data_1_wake','-v7.3')
    end
end

time = data_1_wake.time;
dt = 0.01; % Time resolution of simulation, set in each sim file
downsample_factor = 10; % Set in each sim file
time_index_begin = 10000; % In milliseconds
time_index_end =   11000; % In milliseconds

time_index_range = round(time_index_begin/(dt*downsample_factor)):round(time_index_end/(dt*downsample_factor));

f101 = figure('Color', 'w', 'Units', 'normalized', 'InnerPosition', [0, 0, 1, 1]);
ax1=subplot(5,1,1); plot(time(time_index_range),data_1_wake.PYdr_v(time_index_range,1),'LineWidth',1.5); ylim([-100 50]); set(ax1,'XTickLabel',{},'YTickLabel',{},'box','off')
ax2=subplot(5,1,2); plot(time(time_index_range),data_1_wake.PYso_v(time_index_range,1),'LineWidth',1.5); ylim([-100 50]); set(ax2,'XTickLabel',{},'YTickLabel',{},'box','off')
ax3=subplot(5,1,3); plot(time(time_index_range),data_1_wake.IN_v(time_index_range,1),  'LineWidth',1.5); ylim([-100 50]); set(ax3,'XTickLabel',{},'YTickLabel',{},'box','off')
ax4=subplot(5,1,4); plot(time(time_index_range),data_1_wake.TC_v(time_index_range,1),  'LineWidth',1.5); ylim([-100 50]); set(ax4,'XTickLabel',{},'YTickLabel',{},'box','off')
ax5=subplot(5,1,5); plot(time(time_index_range),data_1_wake.TRN_v(time_index_range,1), 'LineWidth',1.5); ylim([-100 50]); set(ax5,'XTickLabel',{},'YTickLabel',{},'box','off')
print(f101, 'f101_wake_traces', '-dpng','-r300')

% % This plot in particular may require 32GB of RAM due to the amount of
% % spikes
% f102 = dsPlot(data_1_wake, 'plot_type', 'rastergram',...
%               'disregard_analysis_calc',1);
% print(f102, 'f102_wake_raster', '-dpng','-r300')

% % If want awake sim power
% f103 = dsPlot(data_1_wake, 'plot_type', 'power',...
%     'disregard_analysis_calc',1,...
%     'xlim',[0 80]);
% print(f103, 'f103_wake_power', '-dpng')

% Remove the dataset from RAM to free up space
clear data_1_wake

%% Part 2: Simulating, analyzing, and plotting trough-max propofol case

% Load the data if isn't present, and if the data file doesn't exist, run
% the simulation and save its data.
if ~exist('data_2_tmax','var')
    if isfile('data_2_tmax.mat')
        load(strcat(output_dir,'data_2_tmax.mat'))
    else
        % All code and configuration necessary to run and save both the
        % simulation and analyses is inside the following function file.
        % This will produce a folder containing around 17GB of data on
        % disk.
        data_2_tmax = sim_2_tmax_run(output_dir);
        % We'll clear our memory so we can combine and clean all the
        % saved-to-disk datasets
        clear data_2_tmax

        % This will combine and clean the different data sets into a single
        % one, which we will save locally. After this finishes, you can
        % feel free to delete the 'sim_2_tmax_run' folder containing
        % approximately 17GB of data sets. Note that this function
        % requires at LEAST 16GB of RAM available, so it's recommended to
        % close your other programs!
        data_2_tmax = combine_data_files(output_dir, 'sim_2_tmax');
        save(strcat(output_dir,'data_2_tmax.mat'), 'data_2_tmax','-v7.3')      
    end
end

% Choose the portion of the all the time series to display to better
% illustrate the simulation
time = data_2_tmax.time;
dt = 0.01; % Time resolution of simulation, set in each sim file
downsample_factor = 10; % Set in each sim file
time_index_begin = 10000; % In milliseconds
time_index_end =   12000; % In milliseconds
neuron = 1;

time_index_range = round(time_index_begin/(dt*downsample_factor)):round(time_index_end/(dt*downsample_factor));

% Plot the basic cell voltage traces, rastergram, and power spectra
f201 = figure('Color', 'w', 'Units', 'normalized', 'InnerPosition', [0, 0, 1, 1]);
ax1=subplot(5,1,1); plot(time(time_index_range),data_2_tmax.PYdr_v(time_index_range,neuron),'LineWidth',1.5); ylim([-100 50]); set(ax1,'XTickLabel',{},'YTickLabel',{},'box','off')
ax2=subplot(5,1,2); plot(time(time_index_range),data_2_tmax.PYso_v(time_index_range,neuron),'LineWidth',1.5); ylim([-100 50]); set(ax2,'XTickLabel',{},'YTickLabel',{},'box','off')
ax3=subplot(5,1,3); plot(time(time_index_range),data_2_tmax.IN_v(time_index_range,neuron),  'LineWidth',1.5); ylim([-100 50]); set(ax3,'XTickLabel',{},'YTickLabel',{},'box','off')
ax4=subplot(5,1,4); plot(time(time_index_range),data_2_tmax.TC_v(time_index_range,neuron),  'LineWidth',1.5); ylim([-100 50]); set(ax4,'XTickLabel',{},'YTickLabel',{},'box','off')
ax5=subplot(5,1,5); plot(time(time_index_range),data_2_tmax.TRN_v(time_index_range,neuron), 'LineWidth',1.5); ylim([-100 50]); set(ax5,'XTickLabel',{},'YTickLabel',{},'box','off')
print(f201, 'f201_tmax_traces', '-dpng','-r300')

f202 = dsPlot(data_2_tmax, 'plot_type', 'rastergram',...
    'disregard_analysis_calc',1);
print(f202, 'f202_tmax_raster', '-dpng','-r300')

f203 = dsPlot(data_2_tmax, 'plot_type', 'power','xlim',[0 80],...
        'disregard_analysis_calc',1);
print(f203, 'f203_tmax_power', '-dpng','-r300')

% Plot the tmax TCPY and PYPY relevant traces, mean IAMPA currents, and 
% those currents' spectra (our EEG signals)
f204 = figure('Color', 'w', 'Units', 'normalized', 'InnerPosition', [0, 0, 1, 1]);
ax1=subplot(5,2,1);
plot(time(time_index_range),data_2_tmax.PYdr_v(time_index_range,neuron),'LineWidth',1.5); ylim([-100 50]); set(ax1,'XTickLabel',{},'YTickLabel',{},'box','off')
ax2=subplot(5,2,2);
plot(time(time_index_range),data_2_tmax.PYdr_v(time_index_range,neuron),'LineWidth',1.5); ylim([-100 50]); set(ax2,'XTickLabel',{},'YTickLabel',{},'box','off')
ax3=subplot(5,2,3);
plot(time(time_index_range),data_2_tmax.TC_v(time_index_range,neuron),'LineWidth',1.5); ylim([-100 50]); set(ax3,'XTickLabel',{},'YTickLabel',{},'box','off')
% NEIGHBORing PYso
ax4=subplot(5,2,4);
plot(time(time_index_range),data_2_tmax.PYso_v(time_index_range,neuron+1),'LineWidth',1.5); ylim([-100 50]); set(ax4,'XTickLabel',{},'YTickLabel',{},'box','off')
ax5=subplot(5,2,5);
plot(time(time_index_range),mean(data_2_tmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC(time_index_range,:),2),'LineWidth',1.5); set(ax5,'XTickLabel',{},'box','off')
ax6=subplot(5,2,6);
plot(time(time_index_range),mean(data_2_tmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12(time_index_range,:),2),'LineWidth',1.5); set(ax6,'XTickLabel',{},'box','off')
ax7=subplot(5,2,7);
plot(time(time_index_range),mean(data_2_tmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC(time_index_range,:),2),'LineWidth',1.5); set(ax7,'XTickLabel',{},'box','off')
ylim([1.4, 1.6])

ax9=subplot(5,2,9);
plot(data_2_tmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_Power_MUA.frequency(1:1050),pow2db(data_2_tmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_Power_MUA.Pxx(1:1050)),'LineWidth',1.5); xlim([0 80])
ax10=subplot(5,2,10);
plot(data_2_tmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_Power_MUA.frequency(1:1050),pow2db(data_2_tmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_Power_MUA.Pxx(1:1050)),'LineWidth',1.5); xlim([0 80])
print(f204, 'f204_tmax_iampa', '-dpng','-r300')

% Plot the tmax TCPY comodulograms
%
% The filters used for the comodulograms are centered around the Cartesian
% product of {9, 12} Hz for alpha and {0.5, 1.5, 2.5} Hz for SWO.
current_data = data_2_tmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_Comodulograms_MUA.comodulograms;
% 'smaller_data' is for using only a subset of the original comodulogram
% data, since the initial transient of the simulation (the first few
% seconds) takes time to adapt to the steady-state of the simulation.
smaller_data = cell(size(current_data));

f205 = figure('Color', 'w', 'Units', 'normalized', 'InnerPosition', [0, 0, 1, 1]);
for ii=1:size(current_data,1)
    for jj=1:size(current_data,2)
        linear_index = (ii-1)*3+jj;
        smaller_data{ii,jj} = current_data{ii,jj}(:,2:7);
        ax{linear_index} = subplot(2,3,linear_index);
        imagesc(smaller_data{ii,jj})
        set(ax{linear_index}, 'XTickLabel',{},'YTickLabel',{})
        colorbar
    end
end
print(f205, 'f205_tmax_tcpy_comodulograms', '-dpng','-r300')

% Plot the tmax TCPY coupling map
current_data = data_2_tmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_Coupling_MUA;
current_data_amplitudes = current_data.amplitudes;
current_data_angles = current_data.angles;

% 'smaller_data_*' are for removing the coupling-map data pertaining to
% 0.01 Hz-filtered SWO phase, since the power/coupling analysis there is
% skewed by being too close to 0 Hz frequency, and therefore we need to
% throw it out.
smaller_data_amplitudes = current_data_amplitudes(:,2:13);
smaller_data_angles = current_data_angles(:,2:13);

% Note: Don't bother flipping the y-axis up-to-down, since Matlab refuses
% to allow top-to-bottom go highest-to-lowest instead of the other way
% around, ugh. Tried using 'flipud' but no effect.
f206 = figure('Color', 'w');

imagesc(current_data.ph_freq_axis(2:13,:),...
        current_data.ampl_freq_axis,...
        smaller_data_amplitudes)
colorbar
hold on
quiver(current_data.ph_freq_axis(2:13,:),...
       current_data.ampl_freq_axis,...
       cos(smaller_data_angles),...
       sin(smaller_data_angles),...
       'Color', 'k', 'LineWidth', 1.0,'AutoScaleFactor', 0.25)
hold off
print(f206, 'f206_tmax_tcpy_couplingmap', '-dpng','-r300')

% Plot the tmax PYPY comodulograms
% Apparently Matlab truncated my variable name???
current_data = data_2_tmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_Comodulogra.comodulograms;
smaller_data = cell(size(current_data));

f207 = figure('Color', 'w', 'Units', 'normalized', 'InnerPosition', [0, 0, 1, 1]);
for ii=1:size(current_data,1)
    for jj=1:size(current_data,2)
        linear_index = (ii-1)*3+jj;
        smaller_data{ii,jj} = current_data{ii,jj}(:,2:7);
        ax{linear_index} = subplot(2,3,linear_index);
        imagesc(smaller_data{ii,jj})
        set(ax{linear_index}, 'XTickLabel',{},'YTickLabel',{})
        colorbar
    end
end
print(f207, 'f207_tmax_pypy_comodulograms', '-dpng','-r300')

% Plot the tmax PYPY coupling map
current_data = data_2_tmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_Coupling_MU;
current_data_amplitudes = current_data.amplitudes;
current_data_angles = current_data.angles;

smaller_data_amplitudes = current_data_amplitudes(:,2:13);
smaller_data_angles = current_data_angles(:,2:13);

f208 = figure('Color', 'w');

imagesc(current_data.ph_freq_axis(2:13,:),...
        current_data.ampl_freq_axis,...
        smaller_data_amplitudes)
colorbar
hold on
quiver(current_data.ph_freq_axis(2:13,:),...
       current_data.ampl_freq_axis,...
       cos(smaller_data_angles),...
       sin(smaller_data_angles),...
       'Color', 'k', 'LineWidth', 1.0,'AutoScaleFactor', 0.25)
hold off
print(f208, 'f208_tmax_pypy_couplingmap', '-dpng','-r300')

f209 = figure('Color', 'w');

ax1=subplot(2,2,1);
plot(time, data_2_tmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_FiltSigs.filtered_slow{1,1}, 'b', 'LineWidth',1.5)
hold on
plot(time, data_2_tmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_FiltSigs.filtered_slow{1,2}, 'm', 'LineWidth',1.5)
plot(time, data_2_tmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_FiltSigs.filtered_slow{1,3}, 'g', 'LineWidth',1.5)
hold off
set(ax1,'XTickLabel',{},'box','off')
ylim([0 600])

ax2=subplot(2,2,2);
plot(time, data_2_tmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_FiltSigs.filtered_slow{1,1}, 'b', 'LineWidth',1.5)
hold on
plot(time, data_2_tmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_FiltSigs.filtered_slow{1,2}, 'm', 'LineWidth',1.5)
plot(time, data_2_tmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_FiltSigs.filtered_slow{1,3}, 'g', 'LineWidth',1.5)
hold off
set(ax2,'XTickLabel',{},'box','off')
ylim([0 60])

ax3=subplot(2,2,3);
plot(time, data_2_tmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_FiltSigs.filtered_fast{1,1}, 'k', 'LineWidth',1.5)
hold on
plot(time, data_2_tmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_FiltSigs.filtered_fast{1,2}, 'r', 'LineWidth',1.5)
hold off
set(ax3,'XTickLabel',{},'box','off')
ylim([0 600])

ax4=subplot(2,2,4);
plot(time, data_2_tmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_FiltSigs.filtered_fast{1,1}, 'k', 'LineWidth',1.5)
hold on
plot(time, data_2_tmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_FiltSigs.filtered_fast{1,2}, 'r', 'LineWidth',1.5)
hold off
set(ax4,'XTickLabel',{},'box','off')
ylim([0 50])

print(f209, 'f209_tmax_filtsigs', '-dpng','-r300')

f210 = figure('Color', 'w');
ax1=subplot(1,1,1); plot(time(time_index_range),data_2_tmax.TC_PYso_iAMPA_TC_PYso_iAMPA_TC_PYso(time_index_range,neuron),'LineWidth',1.5); set(ax1,'box','off')
print(f210, 'f210_cortthal_current_trace', '-dpng','-r300')

% Remove the dataset from RAM to free up space
clear data_2_tmax

%% Part 3: Simulating, analyzing, and plotting peak-max propofol case

% Load the data if isn't present, and if the data file doesn't exist, run
% the simulation and save its data.
if ~exist('data_3_pmax','var')
    if isfile('data_3_pmax.mat')
        load(strcat(output_dir,'data_3_pmax.mat'))
    else
        % All code and configuration necessary to run and save both the
        % simulation and analyses is inside the following function file.
        % This will produce a folder containing around 17GB of data on
        % disk.
        data_3_pmax = sim_3_pmax_run(output_dir);
        % We'll clear our memory so we can combine and clean all the
        % saved-to-disk datasets
        clear data_3_pmax

        % This will combine and clean the different data sets into a single
        % one, which we will save locally. After this finishes, you can
        % feel free to delete the 'sim_3_pmax_run' folder containing
        % approximately 17GB of data sets. Note that this function
        % requires at LEAST 16GB of RAM available, so it's recommended to
        % close your other programs!
        data_3_pmax = combine_data_files(output_dir, 'sim_3_pmax');
        save(strcat(output_dir,'data_3_pmax.mat'), 'data_3_pmax','-v7.3')      
    end
end

% Choose the portion of the all the time series to display to better
% illustrate the simulation
time = data_3_pmax.time;
dt = 0.01; % Time resolution of simulation, set in each sim file
downsample_factor = 10; % Set in each sim file
time_index_begin =  11001; % In milliseconds
time_index_end =    14000; % In milliseconds
neuron = 2;

time_index_range = round(time_index_begin/(dt*downsample_factor)):round(time_index_end/(dt*downsample_factor));

% Plot the basic cell voltage traces, rastergram, and power spectra
f301 = figure('Color', 'w', 'Units', 'normalized', 'InnerPosition', [0, 0, 1, 1]);
ax1=subplot(5,1,1); plot(time(time_index_range),data_3_pmax.PYdr_v(time_index_range,neuron),'LineWidth',1.5); ylim([-100 50]); set(ax1,'XTickLabel',{},'YTickLabel',{},'box','off')
ax2=subplot(5,1,2); plot(time(time_index_range),data_3_pmax.PYso_v(time_index_range,neuron),'LineWidth',1.5); ylim([-100 50]); set(ax2,'XTickLabel',{},'YTickLabel',{},'box','off')
ax3=subplot(5,1,3); plot(time(time_index_range),data_3_pmax.IN_v(time_index_range,neuron),  'LineWidth',1.5); ylim([-100 50]); set(ax3,'XTickLabel',{},'YTickLabel',{},'box','off')
ax4=subplot(5,1,4); plot(time(time_index_range),data_3_pmax.TC_v(time_index_range,neuron),  'LineWidth',1.5); ylim([-100 50]); set(ax4,'XTickLabel',{},'YTickLabel',{},'box','off')
ax5=subplot(5,1,5); plot(time(time_index_range),data_3_pmax.TRN_v(time_index_range,neuron), 'LineWidth',1.5); ylim([-100 50]); set(ax5,'XTickLabel',{},'YTickLabel',{},'box','off')
print(f301, 'f301_pmax_traces', '-dpng','-r300')

f302 = dsPlot(data_3_pmax, 'plot_type', 'rastergram',...
    'disregard_analysis_calc',1);
print(f302, 'f302_pmax_raster', '-dpng','-r300')

f303 = dsPlot(data_3_pmax, 'plot_type', 'power','xlim',[0 80],...
        'disregard_analysis_calc',1);
print(f303, 'f303_pmax_power', '-dpng','-r300')

% Plot the pmax TCPY and PYPY relevant traces, mean IAMPA currents, and 
% those currents' spectra (our EEG signals)
f304 = figure('Color', 'w', 'Units', 'normalized', 'InnerPosition', [0, 0, 1, 1]);
ax1=subplot(4,2,1);
plot(time(time_index_range),data_3_pmax.PYdr_v(time_index_range,neuron),'LineWidth',1.5); ylim([-100 50]); set(ax1,'XTickLabel',{},'YTickLabel',{},'box','off')
ax2=subplot(4,2,2);
plot(time(time_index_range),data_3_pmax.PYdr_v(time_index_range,neuron),'LineWidth',1.5); ylim([-100 50]); set(ax2,'XTickLabel',{},'YTickLabel',{},'box','off')
ax3=subplot(4,2,3);
plot(time(time_index_range),data_3_pmax.TC_v(time_index_range,neuron),  'LineWidth',1.5); ylim([-100 50]); set(ax3,'XTickLabel',{},'YTickLabel',{},'box','off')
% NEIGHBORing PYso
ax4=subplot(4,2,4);
plot(time(time_index_range),data_3_pmax.PYso_v(time_index_range,neuron+1),'LineWidth',1.5); ylim([-100 50]); set(ax4,'XTickLabel',{},'YTickLabel',{},'box','off')
ax5=subplot(4,2,5);
plot(time(time_index_range),mean(data_3_pmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC(time_index_range,:),2),'LineWidth',1.5); set(ax5,'XTickLabel',{},'box','off')
ax6=subplot(4,2,6);
plot(time(time_index_range),mean(data_3_pmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12(time_index_range,:),2),'LineWidth',1.5); set(ax6,'XTickLabel',{},'box','off')
ax7=subplot(4,2,7);
plot(data_3_pmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_Power_MUA.frequency(1:1050),pow2db(data_3_pmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_Power_MUA.Pxx(1:1050)),'LineWidth',1.5); xlim([0 80])
ax8=subplot(4,2,8);
plot(data_3_pmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_Power_MUA.frequency(1:1050),pow2db(data_3_pmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_Power_MUA.Pxx(1:1050)),'LineWidth',1.5); xlim([0 80])
print(f304, 'f304_pmax_iampa', '-dpng','-r300')

% Plot the pmax TCPY comodulograms
%
% The filters used for the comodulograms are centered around the Cartesian
% product of {9, 12} Hz for alpha and {0.01, 1.01, 2.01} Hz for SWO.
current_data = data_3_pmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_Comodulograms_MUA.comodulograms;
% 'smaller_data' is for using only a subset of the original comodulogram
% data, since the initial transient of the simulation (the first few
% seconds) takes time to adapt to the steady-state of the simulation.
smaller_data = cell(size(current_data));

f305 = figure('Color', 'w', 'Units', 'normalized', 'InnerPosition', [0, 0, 1, 1]);
for ii=1:size(current_data,1)
    for jj=1:size(current_data,2)
        linear_index = (ii-1)*3+jj;
        smaller_data{ii,jj} = current_data{ii,jj}(:,2:7);
        ax{linear_index} = subplot(2,3,linear_index);
        imagesc(smaller_data{ii,jj})
        set(ax{linear_index}, 'XTickLabel',{},'YTickLabel',{})
        colorbar
    end
end
print(f305, 'f305_pmax_tcpy_comodulograms', '-dpng','-r300')

% Plot the pmax TCPY coupling map
current_data = data_3_pmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_Coupling_MUA;
current_data_amplitudes = current_data.amplitudes;
current_data_angles = current_data.angles;

% 'smaller_data_*' are for removing the coupling-map data pertaining to
% 0.01 Hz-filtered SWO phase, since the power/coupling analysis there is
% skewed by being too close to 0 Hz frequency, and therefore we need to
% throw it out.
smaller_data_amplitudes = current_data_amplitudes(:,2:13);
smaller_data_angles = current_data_angles(:,2:13);

% Note: Don't bother flipping the y-axis up-to-down, since Matlab refuses
% to allow top-to-bottom go highest-to-lowest instead of the other way
% around, ugh. Tried using 'flipud' but no effect.
f306 = figure('Color', 'w');

imagesc(current_data.ph_freq_axis(2:13,:),...
        current_data.ampl_freq_axis,...
        smaller_data_amplitudes)
colorbar
hold on
quiver(current_data.ph_freq_axis(2:13,:),...
       current_data.ampl_freq_axis,...
       cos(smaller_data_angles),...
       sin(smaller_data_angles),...
       'Color', 'k', 'LineWidth', 1.0,'AutoScaleFactor', 0.25)
hold off
print(f306, 'f306_pmax_tcpy_couplingmap', '-dpng','-r300')

% Plot the pmax PYPY comodulograms
% Apparently Matlab truncated my variable name???
current_data = data_3_pmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_Comodulogra.comodulograms;
smaller_data = cell(size(current_data));

f307 = figure('Color', 'w', 'Units', 'normalized', 'InnerPosition', [0, 0, 1, 1]);
for ii=1:size(current_data,1)
    for jj=1:size(current_data,2)
        linear_index = (ii-1)*3+jj;
        smaller_data{ii,jj} = current_data{ii,jj}(:,2:7);
        ax{linear_index} = subplot(2,3,linear_index);
        imagesc(smaller_data{ii,jj})
        set(ax{linear_index}, 'XTickLabel',{},'YTickLabel',{})
        colorbar
    end
end
print(f307, 'f307_pmax_pypy_comodulograms', '-dpng','-r300')

% Plot the pmax PYPY coupling map
current_data = data_3_pmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_Coupling_MU;
current_data_amplitudes = current_data.amplitudes;
current_data_angles = current_data.angles;

smaller_data_amplitudes = current_data_amplitudes(:,2:13);
smaller_data_angles = current_data_angles(:,2:13);

f308 = figure('Color', 'w');

imagesc(current_data.ph_freq_axis(2:13,:),...
        current_data.ampl_freq_axis,...
        smaller_data_amplitudes)
colorbar
hold on
quiver(current_data.ph_freq_axis(2:13,:),...
       current_data.ampl_freq_axis,...
       cos(smaller_data_angles),...
       sin(smaller_data_angles),...
       'Color', 'k', 'LineWidth', 1.0,'AutoScaleFactor', 0.25)
hold off
print(f308, 'f308_pmax_pypy_couplingmap', '-dpng','-r300')

f309 = figure('Color', 'w');

ax1=subplot(2,2,1);
plot(time, data_3_pmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_FiltSigs.filtered_slow{1,1}, 'b', 'LineWidth',1.5)
hold on
plot(time, data_3_pmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_FiltSigs.filtered_slow{1,2}, 'm', 'LineWidth',1.5)
plot(time, data_3_pmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_FiltSigs.filtered_slow{1,3}, 'g', 'LineWidth',1.5)
hold off
set(ax1,'XTickLabel',{},'box','off')
ylim([0 600])

ax2=subplot(2,2,2);
plot(time, data_3_pmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_FiltSigs.filtered_slow{1,1}, 'b', 'LineWidth',1.5)
hold on
plot(time, data_3_pmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_FiltSigs.filtered_slow{1,2}, 'm', 'LineWidth',1.5)
plot(time, data_3_pmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_FiltSigs.filtered_slow{1,3}, 'g', 'LineWidth',1.5)
hold off
set(ax2,'XTickLabel',{},'box','off')
ylim([0 60])

ax3=subplot(2,2,3);
plot(time, data_3_pmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_FiltSigs.filtered_fast{1,1}, 'k', 'LineWidth',1.5)
hold on
plot(time, data_3_pmax.PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC_FiltSigs.filtered_fast{1,2}, 'r', 'LineWidth',1.5)
hold off
set(ax3,'XTickLabel',{},'box','off')
ylim([0 600])

ax4=subplot(2,2,4);
plot(time, data_3_pmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_FiltSigs.filtered_fast{1,1}, 'k', 'LineWidth',1.5)
hold on
plot(time, data_3_pmax.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_FiltSigs.filtered_fast{1,2}, 'r', 'LineWidth',1.5)
hold off
set(ax4,'XTickLabel',{},'box','off')
ylim([0 50])

print(f309, 'f309_pmax_filtsigs', '-dpng','-r300')

f310 = figure('Color', 'w');
ax1=subplot(1,1,1); plot(time(time_index_range),data_3_pmax.TC_PYso_iAMPA_TC_PYso_iAMPA_TC_PYso(time_index_range,neuron),'LineWidth',1.5); set(ax1,'box','off')
print(f310, 'f310_cortthal_current_trace', '-dpng','-r300')

% Remove the dataset from RAM to free up space
clear data_3_pmax

%% Part 4: Simulating, analyzing, and plotting hi-GABAA, hi-ACh case

if ~exist('data_4_GA_ACh','var')
    if isfile('data_4_GA_ACh.mat')
        % Try to load the simulation data, if you've already run the sim
        % and it exists.
        load(strcat(output_dir,'data_4_GA_ACh.mat'))
    else
        % All code and configuration necessary to run and save both the
        % simulation and analyses is inside the following function file.
        % This will produce a folder containing around 17GB of data on
        % disk.
        data_4_GA_ACh = sim_4_GA_ACh_run(output_dir);
        % We'll clear our memory so we can combine and clean all the
        % saved-to-disk datasets
        clear data_4_GA_ACh

        % This will combine and clean the different data sets into a single
        % one, which we will save locally. After this finishes, you can
        % feel free to delete the 'sim_1_wake_run' folder containing
        % approximately 17GB of data sets. Note that this function
        % requires at LEAST 16GB of RAM available, so it's recommended to
        % close your other programs!
        data_4_GA_ACh = combine_data_files(output_dir, 'sim_4_GA_ACh');
        save(strcat(output_dir,'data_4_GA_ACh.mat'), 'data_4_GA_ACh','-v7.3')
    end
end

time = data_4_GA_ACh.time;
dt = 0.01; % Time resolution of simulation, set in each sim file
downsample_factor = 10; % Set in each sim file
time_index_begin = 10001; % In milliseconds
time_index_end =   12000; % In milliseconds
neuron = 20;

time_index_range = round(time_index_begin/(dt*downsample_factor)):round(time_index_end/(dt*downsample_factor));

f401 = figure('Color', 'w', 'Units', 'normalized', 'InnerPosition', [0, 0, 1, 1]);
ax1=subplot(5,1,1); plot(time(time_index_range),data_4_GA_ACh.PYdr_v(time_index_range,neuron),'LineWidth',1.5); ylim([-100 50]); set(ax1,'XTickLabel',{},'YTickLabel',{},'box','off')
ax2=subplot(5,1,2); plot(time(time_index_range),data_4_GA_ACh.PYso_v(time_index_range,neuron),'LineWidth',1.5); ylim([-100 50]); set(ax2,'XTickLabel',{},'YTickLabel',{},'box','off')
ax3=subplot(5,1,3); plot(time(time_index_range),data_4_GA_ACh.IN_v(time_index_range,neuron),  'LineWidth',1.5); ylim([-100 50]); set(ax3,'XTickLabel',{},'YTickLabel',{},'box','off')
ax4=subplot(5,1,4); plot(time(time_index_range),data_4_GA_ACh.TC_v(time_index_range,neuron),  'LineWidth',1.5); ylim([-100 50]); set(ax4,'XTickLabel',{},'YTickLabel',{},'box','off')
ax5=subplot(5,1,5); plot(time(time_index_range),data_4_GA_ACh.TRN_v(time_index_range,neuron), 'LineWidth',1.5); ylim([-100 50]); set(ax5,'XTickLabel',{},'YTickLabel',{},'box','off')
print(f401, 'f401_GA_ACh_traces', '-dpng','-r300')

% Remove the dataset from RAM to free up space
clear data_4_GA_ACh

%% Part 5

if ~exist('data_5_hiGA','var')
    if isfile('data_5_hiGA.mat')
        % Try to load the simulation data, if you've already run the sim
        % and it exists.
        load(strcat(output_dir,'data_5_hiGA.mat'))
    else
        % All code and configuration necessary to run and save both the
        % simulation and analyses is inside the following function file.
        % This will produce a folder containing around 17GB of data on
        % disk.
        data_5_hiGA = sim_5_hiGA_run(output_dir);
        % We'll clear our memory so we can combine and clean all the
        % saved-to-disk datasets
        clear data_5_hiGA

        % This will combine and clean the different data sets into a single
        % one, which we will save locally. After this finishes, you can
        % feel free to delete the 'sim_1_wake_run' folder containing
        % approximately 17GB of data sets. Note that this function
        % requires at LEAST 16GB of RAM available, so it's recommended to
        % close your other programs!
        data_5_hiGA = combine_data_files(output_dir, 'sim_5_hiGA');
        save(strcat(output_dir,'data_5_hiGA.mat'), 'data_5_hiGA','-v7.3')
    end
end

time = data_5_hiGA.time;
dt = 0.01; % Time resolution of simulation, set in each sim file
downsample_factor = 10; % Set in each sim file
time_index_begin = 10001; % In milliseconds
time_index_end =   14000; % In milliseconds
neuron = 5;

time_index_range = round(time_index_begin/(dt*downsample_factor)):round(time_index_end/(dt*downsample_factor));

f501 = figure('Color', 'w', 'Units', 'normalized', 'InnerPosition', [0, 0, 1, 1]);
ax1=subplot(3,1,1); plot(time(time_index_range),data_5_hiGA.PYdr_v(time_index_range,neuron),'LineWidth',1.5); ylim([-100 50]); set(ax1,'XTickLabel',{},'YTickLabel',{},'box','off')
ax2=subplot(3,1,2); plot(time(time_index_range),data_5_hiGA.PYso_v(time_index_range,neuron),'LineWidth',1.5); ylim([-100 50]); set(ax2,'XTickLabel',{},'YTickLabel',{},'box','off')
ax3=subplot(3,1,3); plot(time(time_index_range),data_5_hiGA.IN_v(time_index_range,neuron),  'LineWidth',1.5); ylim([-100 50]); set(ax3,'XTickLabel',{},'YTickLabel',{},'box','off')
print(f501, 'f501_hiGA_traces', '-dpng','-r300')

f502 = dsPlot(data_5_hiGA, 'plot_type', 'rastergram',...
    'disregard_analysis_calc',1);
print(f502, 'f502_hiGA_raster', '-dpng','-r300')

f503 = dsPlot(data_5_hiGA, 'plot_type', 'power','xlim',[0 80],...
        'disregard_analysis_calc',1);
print(f503, 'f505_hiGA_power', '-dpng','-r300')

% Plot the pmax TCPY and PYPY relevant traces, mean IAMPA currents, and 
% those currents' spectra (our EEG signals)
f504 = figure('Color', 'w', 'Units', 'normalized', 'InnerPosition', [0, 0, 1, 1]);
ax1=subplot(4,1,1);
plot(time(time_index_range),data_5_hiGA.PYdr_v(time_index_range,neuron),'LineWidth',1.5); ylim([-100 50]); set(ax1,'XTickLabel',{},'YTickLabel',{},'box','off')
% NEIGHBORing PYso
ax2=subplot(4,1,2);
plot(time(time_index_range),data_5_hiGA.PYso_v(time_index_range,neuron+1),'LineWidth',1.5); ylim([-100 50]); set(ax2,'XTickLabel',{},'YTickLabel',{},'box','off')
ax3=subplot(4,1,3);
plot(time(time_index_range),mean(data_5_hiGA.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12(time_index_range,:),2),'LineWidth',1.5); set(ax3,'XTickLabel',{},'box','off')
ax4=subplot(4,1,4);
plot(data_5_hiGA.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_Power_MUA.frequency(1:1050),pow2db(data_5_hiGA.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_Power_MUA.Pxx(1:1050)),'LineWidth',1.5); xlim([0 80])
print(f504, 'f504_hiGA_iampa', '-dpng','-r300')

% Plot the pmax PYPY comodulograms
% Apparently Matlab truncated my variable name???
current_data = data_5_hiGA.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_Comodulogra.comodulograms;
smaller_data = cell(size(current_data));

f505 = figure('Color', 'w', 'Units', 'normalized', 'InnerPosition', [0, 0, 1, 1]);
for ii=1:size(current_data,1)
    for jj=1:size(current_data,2)
        linear_index = (ii-1)*3+jj;
        smaller_data{ii,jj} = current_data{ii,jj}(:,2:7);
        ax{linear_index} = subplot(2,3,linear_index);
        imagesc(smaller_data{ii,jj})
        set(ax{linear_index}, 'XTickLabel',{},'YTickLabel',{})
        colorbar
    end
end
print(f505, 'f505_hiGA_pypy_comodulograms', '-dpng','-r300')

% Plot the pmax PYPY coupling map
current_data = data_5_hiGA.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_Coupling_MU;
current_data_amplitudes = current_data.amplitudes;
current_data_angles = current_data.angles;

smaller_data_amplitudes = current_data_amplitudes(:,2:13);
smaller_data_angles = current_data_angles(:,2:13);

f506 = figure('Color', 'w');

imagesc(current_data.ph_freq_axis(2:13,:),...
        current_data.ampl_freq_axis,...
        smaller_data_amplitudes)
colorbar
hold on
quiver(current_data.ph_freq_axis(2:13,:),...
       current_data.ampl_freq_axis,...
       cos(smaller_data_angles),...
       sin(smaller_data_angles),...
       'Color', 'k', 'LineWidth', 1.0,'AutoScaleFactor', 0.25)
hold off
print(f506, 'f506_hiGA_pypy_couplingmap', '-dpng','-r300')

% Plot IAMPA filtered signals
f507 = figure('Color', 'w');

ax1=subplot(2,1,1);
plot(time, data_5_hiGA.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_FiltSigs.filtered_slow{1,1}, 'b', 'LineWidth',1.5)
hold on
plot(time, data_5_hiGA.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_FiltSigs.filtered_slow{1,2}, 'm', 'LineWidth',1.5)
plot(time, data_5_hiGA.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_FiltSigs.filtered_slow{1,3}, 'g', 'LineWidth',1.5)
hold off
set(ax1,'XTickLabel',{},'box','off')

ax2=subplot(2,1,2);
plot(time, data_5_hiGA.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_FiltSigs.filtered_fast{1,1}, 'k', 'LineWidth',1.5)
hold on
plot(time, data_5_hiGA.PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12_FiltSigs.filtered_fast{1,2}, 'r', 'LineWidth',1.5)
hold off
set(ax4,'XTickLabel',{},'box','off')

print(f507, 'f507_hiGA_filtsigs', '-dpng','-r300')

% % Remove the dataset from RAM to free up space
% clear data_5_hiGA
