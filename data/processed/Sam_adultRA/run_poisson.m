function run_poisson(case_ids)


for i = 1:length(case_ids)
    case_id = case_ids{i};
    
    
    [q_maxs_poisson_test, H_maxs_poisson_test] = poisson_analysis(case_id);
    
    save(case_id, 'q_maxs_poisson_test','H_maxs_poisson_test','-append');
    
end

end