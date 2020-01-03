#!/bin/bash

# This script runs AutoGrow4 from within the docker container
# It serves as the ENTRYPOINT of the docker.

if echo $* | grep -q "test"; then
  exit;
fi

# Handle moving old run from UserFiles/ to /Outputfolder/ if needed
path_old_runs=UserFiles/old_runs/Run_0
if test -d "$path_old_runs"; then


    # If previous run was performed in docker the generations may be compressed
    # In which case we need to unzip the Outputfolder.zip file
    zip_file=/UserFiles/old_runs/Run_0/Outputfolder.zip
    if test -f "$zip_file"; then
        echo "#########################################"
        echo "BASING AUTOGROW RUN EXCLUSIVELY ON Outputfolder.zip FROM PREVIOUS RUN"
        echo "IF THERE IS INFORMATION NOT CONTAINED IN Outputfolder.zip FROM PREVIOUS RUN ATTEMPTS"
        echo "PLEASE REMOVE THE Outputfolder.zip FROM THE DIRECTORY AND RESTART AUTOGROW RUN"
        echo "#########################################"
        
        mkdir /Outputfolder/temp
        unzip $zip_file -d /Outputfolder/temp/ >> zip_unzip_record.txt

        # mv all previous run data to the Outputfolder/ (this will contain generations info)
        mv Outputfolder/temp/Outputfolder/* Outputfolder/

        # Delete unnecessary files
        rm -rf UserFiles/old_runs
        rm -rf /Outputfolder/temp

    else
        mkdir /Outputfolder/Run_0
        # This previous run does not contain a zip file
        cp -r $path_old_runs/* /Outputfolder/Run_0
    fi
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

    # Add spacer on log files
    if test -f "$output_file"; then
        echo "" >> $output_file
        echo "" >> $output_file
        echo "" >> $output_file
        echo "###################################\n" >> $output_file
        echo "CONTINUE PREVIOUS AUTOGROW4 RUN AT: " >> $output_file
        date >> $output_file
        echo "###################################\n" >> $output_file
        echo "" >> $output_file
        echo "" >> $output_file
        echo "" >> $output_file
    fi
    if test -f "$error_file"; then
        echo "" >> $error_file
        echo "" >> $error_file
        echo "" >> $error_file
        echo "###################################\n" >> $error_file
        echo "CONTINUE PREVIOUS AUTOGROW4 RUN AT: " >> $error_file
        date >> $error_file
        echo "###################################\n" >> $error_file
        echo "" >> $error_file
        echo "" >> $error_file
        echo "" >> $error_file
    fi

    echo "/root/miniconda3/bin/python autogrow4/RunAutogrow.py -j /UserFiles/docker_json_vars.json >> $output_file 2> $error_file"
    # Run autogrow
    /root/miniconda3/bin/python autogrow4/RunAutogrow.py \
        -j /UserFiles/docker_json_vars.json >> $output_file 2>> $error_file

    echo "Completed AutoGrow4 Run"
    # Zip up results and move to working directory.
    # cd //Output/
    chmod -R a+rwx Outputfolder/
    zip -r Outputfolder.zip Outputfolder/ >> zip_unzip_record.txt
    chmod -R a+rwx Outputfolder.zip
fi
# For interactive
bash
