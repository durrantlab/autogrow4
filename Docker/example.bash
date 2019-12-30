#!/bin/bash

SCRIPT_PATH="$(readlink -f "$0")"
DIR_PATH=`dirname "$SCRIPT_PATH"`

# Make sure we are in the autogrow4/Docker/ directory
cd $DIR_PATH

# sudo should only be run in Linux or MacOS
# If Windows please instead just open the terminal with admin priveledges
#   and omit the 'sudo'

sudo python ./autogrow_in_docker.py -j ./sample_autogrow_json.json
