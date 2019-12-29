#!/bin/bash

SCRIPT_PATH="$(readlink -f "$0")"
DIR_PATH=`dirname "$SCRIPT_PATH"`
# Change into top of the autogrow directory
# cd $DIR_PATH
# cd ../
if echo $@ | grep -q "test"; then
  exit;
fi

if echo $1 | grep -q "json"; then

    function my_date {
    date "+%y_%m_%d"
    }


    date_time=$(my_date)
    error=_error.txt
    output=_output.txt
    Outputfolder=/Outputfolder/
    output_file=$Outputfolder$date_time$output
    error_file=$Outputfolder$date_time$error
    echo "/root/miniconda3/bin/python autogrow4/RunAutogrow.py -j $1 >> $output_file 2> $error_file"
    # Run autogrow
    /root/miniconda3/bin/python autogrow4/RunAutogrow.py \
        -j $1 >> $output_file 2> $error_file

fi

# Zip up results and move to working directory.
# cd /autogrow_work_dir/
# zip -r autogrow_output.zip autogrow_output/

# For interactive
bash