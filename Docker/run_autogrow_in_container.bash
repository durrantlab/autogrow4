#!/bin/bash

SCRIPT_PATH="$(readlink -f "$0")"
DIR_PATH=`dirname "$SCRIPT_PATH"`
# Change into top of the autogrow directory
# cd $DIR_PATH
# cd ../
if echo $* | grep -q "test"; then
  exit;
fi
if echo $* | grep -q "json"; then
    function my_date {
    date "+%y_%m_%d"
    }

    echo Running AutoGrow4

    date_time=$(my_date)
    error=_error.txt
    output=_output.txt
    Outputfolder=/Outputfolder/
    output_file=$Outputfolder$date_time$output
    error_file=$Outputfolder$date_time$error
    echo "/root/miniconda3/bin/python autogrow4/RunAutogrow.py -j /UserFiles/docker_json_vars.json >> $output_file 2> $error_file"
    # Run autogrow
    /root/miniconda3/bin/python autogrow4/RunAutogrow.py \
        -j /UserFiles/docker_json_vars.json >> $output_file 2> $error_file

    # Zip up results and move to working directory.
    # cd //Output/
    chmod -R a+rwx Outputfolder/
    zip -r Outputfolder.zip Outputfolder/
    chmod -R a+rwx Outputfolder.zip
fi
# For interactive
bash