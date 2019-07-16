function combined_data = combine_data_files(output_dir, sim_name)
%COMBINE_DATA_FILES This function

combined_data = getfield(load(strcat(output_dir, sim_name,...
    '_run/data/study_sim1_analysis1_dsCalcFR.mat')), 'result');

if strcmp(sim_name, 'sim_5_hiGA')
    other_data_names = {...
        '2_dsCalcPower',...
        '3_dsCalcPower', '4_dsCalcPACFiltSignals',...
        '5_dsCalcComodulograms','6_dsCalcCoupling',...
        };
else
    other_data_names = {...
        '2_dsCalcPower',...
        '3_dsCalcPower', '4_dsCalcPACFiltSignals',...
        '5_dsCalcComodulograms','6_dsCalcCoupling',...
        '7_dsCalcPower', '8_dsCalcPACFiltSignals',...
        '9_dsCalcComodulograms','10_dsCalcCoupling',...
        };
end

for ii=1:length(other_data_names)
    other_data = getfield(load(strcat(output_dir, sim_name,'_run/data/study_sim1_analysis',...
        other_data_names{ii}, '.mat')), 'result');
    
    fields_combined = fieldnames(combined_data);
    fields_other = fieldnames(other_data);
    
    % Find what fields the structs have in common
    common_fields = intersect(fields_combined, fields_other);
    
    % Remove the identical fields from the secondary data struct
    for jj=1:length(common_fields)
        other_data = rmfield(other_data, common_fields{jj});
    end
    
    % Now add the unique fields in the secondary data struct to our
    % combined struct
    fields_other_pruned = fieldnames(other_data);
    for kk=1:length(fields_other_pruned)
        combined_data.(fields_other_pruned{kk}) = other_data.(fields_other_pruned{kk});
    end

    clear other_data
end

end
