#!/bin/bash
temp_user_files=$1
outfolder=$2
echo "###################################"
echo $temp_user_files
echo $outfolder
echo "###################################"
# sudo docker run --rm  -v $temp_user_files \
#     --name autogrow
sudo docker run autogrow4 \
    --name autogrow4 -/UserFiles/docker_json_vars.json    


CONTAINER_ID=$(sudo docker ps -alq)
echo "###################################"
echo "In THIS CONTAINER_ID: "$CONTAINER_ID
echo "###################################"
sudo docker cp $CONTAINER_ID:/Outputfolder.zip $outfolder
