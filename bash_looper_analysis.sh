input="target_list.txt"
while IFS= read -r line
do
  #echo $line"_prior.txt" > run_out_puts/$line"_prior.txt" &
  python3 -u main_analysis.py -s $line -p False -n 1 -c Ruprecht_147  &> run_out_puts/$line"no_prior.txt" &
done < "$input"
