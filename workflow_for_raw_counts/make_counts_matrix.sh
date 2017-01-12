echo "path to sample: $1"
echo "path to output_file: $2"


gzip -d $1*gz
python shift_atac_tagalign.py $1*tagAlign $1/slopped
bedtools coverage -counts -a ppr.merged.bed -b $1/slopped | cut -f4 > $2
