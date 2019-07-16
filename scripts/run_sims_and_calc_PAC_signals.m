%{
run_sims_and_calc_PAC_signals - If you don't have appropriately-named data
files located in the current directory (e.g. you're running the simulations
for the first time), then this runs each simulation+analysis script,
combines the relevant analyses data structures, then saves it to disk,
then in certain cases (sims 1-3) runs an additional analysis before saving
the final version of the data file to disk.

This will obviously take a long time to run, see the README for details.

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

data_1_wake = dsCalcPACFiltSignals(data_1_wake, 'variable', 'PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12');
data_1_wake = dsCalcPACFiltSignals(data_1_wake, 'variable', 'PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC');

save(strcat(output_dir,'data_1_wake.mat'), 'data_1_wake','-v7.3')

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

data_2_tmax = dsCalcPACFiltSignals(data_2_tmax, 'variable', 'PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12');
data_2_tmax = dsCalcPACFiltSignals(data_2_tmax, 'variable', 'PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC');

save(strcat(output_dir,'data_2_tmax.mat'), 'data_2_tmax','-v7.3')

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

data_3_pmax = dsCalcPACFiltSignals(data_3_pmax, 'variable', 'PYdr_PYso_iAMPA_PYdr_PYso_JB12_IAMPA_PYdr_PYso_JB12');
data_3_pmax = dsCalcPACFiltSignals(data_3_pmax, 'variable', 'PYdr_TC_iAMPA_PYdr_TC_iAMPA_PYdr_TC');

save(strcat(output_dir,'data_3_pmax.mat'), 'data_3_pmax','-v7.3')

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

% Remove the dataset from RAM to free up space
clear data_5_hiGA