function run_poisson2(case_ids)


for i = 1:length(case_ids)
    case_id = case_ids{i};
    
    
    [q_maxs_poisson_test2, H_maxs_poisson_test2] = poisson_analysis2(case_id);
    
    save(case_id, 'q_maxs_poisson_test2','H_maxs_poisson_test2','-append');
    
end

end