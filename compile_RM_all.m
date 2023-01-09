function compile_RM_all(num_subj_per_batch, num_batch)
for version = [1:9,11,13,73]
    compile_RM_nonparam(version, num_subj_per_batch, num_batch);
end

compile_RM_param(1, true, num_subj_per_batch, num_batch); compile_RM_param(1, false, num_subj_per_batch, num_batch);
compile_RM_param(3, true, num_subj_per_batch, num_batch); compile_RM_param(3, false, num_subj_per_batch, num_batch);
compile_RM_param(4, true, num_subj_per_batch, num_batch);
compile_RM_param(5, true, num_subj_per_batch, num_batch); compile_RM_param(5, false, num_subj_per_batch, num_batch);
compile_RM_param(7, true, num_subj_per_batch, num_batch); compile_RM_param(7, false, num_subj_per_batch, num_batch);
compile_RM_param(9, true, num_subj_per_batch, num_batch); compile_RM_param(9, false, num_subj_per_batch, num_batch);
compile_RM_param(11, true, num_subj_per_batch, num_batch);
compile_RM_param(73, true, num_subj_per_batch, num_batch);
end