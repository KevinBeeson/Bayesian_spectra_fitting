OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1 
MKL_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1
NUMEXPR_NUM_THREADS=1
input="target_list.txt"
while IFS= read -r line
do
  #echo $line"_prior.txt" > run_out_puts/$line"_prior.txt" &
  python3 -u main_analysis.py -s $line -p False -n 2 -c Ruprecht_147  &> run_out_puts/$line"_no_normalization_no_prior.txt" &
done < "$input"
