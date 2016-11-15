sh make_counts_matrix.sh earlyG1_R1_controls/out/align/pooled_pseudo_reps/ppr1/ sample1
sh make_counts_matrix.sh earlyG1_R1_controls/out/align/pooled_pseudo_reps/ppr2/ sample2
sh make_counts_matrix.sh earlyG1_R1_treated/out/align/pooled_pseudo_reps/ppr1/ sample3
sh make_counts_matrix.sh earlyG1_R1_treated/out/align/pooled_pseudo_reps/ppr2/ sample4
paste sample1 sample2 sample3 sample4 > counts.txt

