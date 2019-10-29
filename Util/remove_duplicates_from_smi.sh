
# filename=/home/jacob/Downloads/zinc15_available/pass_sanitize.smi
# awk '!seen[$1]++' $filename >> output.txt
# awk '!seen[$2]++' output.txt >> output2.txt

filename=/home/jacob/Documents/autogrow4/source_compounds/ZINC_fragments.smi
tmp_str=_tmp_dummy.smi
no_dup_str=_no_dup.smi
tmp_filename=$filename$tmp_str
new_filename=$filename$no_dup_str
echo $tmp_filename

awk '!seen[$1]++' $filename >> $tmp_filename
echo $new_filename
awk '!seen[$2]++' $tmp_filename >> $new_filename

rm $tmp_filename
