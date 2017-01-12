sh make_counts_matrix.sh treated/out/align/pooled_pseudo_reps/ppr1/ treated1
sh make_counts_matrix.sh treated/out/align/pooled_pseudo_reps/ppr2/ treated2
sh make_counts_matrix.sh control/out/align/pooled_pseudo_reps/ppr1/ control1
sh make_counts_matrix.sh control/out/align/pooled_pseudo_reps/ppr2/ control2
paste control1 control2 treated1 treated2 > counts.txt


