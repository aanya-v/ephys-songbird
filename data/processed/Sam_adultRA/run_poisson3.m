function run_poisson3(case_ids)


for i = 1:length(case_ids)
    case_id = case_ids{i};
    
    
    [q_maxs_poisson_test3, H_maxs_poisson_test3] = poisson_analysis2(case_id);
    
    save(case_id, 'q_maxs_poisson_test3','H_maxs_poisson_test3','-append');
    
end

end