function bayesinfer_output_files2mat(out_file, bene_score, path_isswoh, path_bsswh, path_iswhc, path_iswohc,path_bswhc, path_stats)
% input:
% out_file   : string containig result file name,result_results_n_p_a_s.mat
% path_isswoh : path to file initial_state_scores_wihtout_hidden
% path_bsswh : path to file best_state_scores_with_hidden
% path_iswhc : path to file initial_state_with_hidden_counts
% path_iswohc: path to file initial_state_without_hidden_counts
% path_bswhc : path to file best_state_with_hidden_counts
% path_stats : path to file statistics.mat
% 

isswoh   =dlmread(path_isswoh);
bsswh   =dlmread(path_bsswh);
iswhc   =dlmread(path_iswhc);
iswohc  =dlmread(path_iswohc);
bswhc   =dlmread(path_bswhc);
load(path_stats);

save(out_file, 'bene_score', 'isswoh', 'bsswh', 'iswhc', 'iswohc','bswhc', 'statistics', 'hidden_found_fci', 'pdag_tmp', 'pdag_with_hidden', 'bnet');
