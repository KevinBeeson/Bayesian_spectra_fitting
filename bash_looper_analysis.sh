input="melotte_22_targets.txt"
while IFS= read -r line
do
  #echo $line"_prior.txt" > run_out_puts/$line"_prior.txt" &
  python3 main_analysis.py -s $line -p False -n 2 > run_out_puts/$line"_prior.txt" &
done < "$input"
