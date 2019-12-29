# #!/bin/bash

# dir_name=$0
# echo $dir_name
# SCRIPT_PATH="$(readlink -f "$0")"
# echo $SCRIPT_PATH

# # dir_name "$(readlink -f "$0")"
# # echo $dir_name
# echo ""
# dir_name=`dirname "$SCRIPT_PATH"`
# echo $dir_name


function my_date {
  date "+%y_%m_%d"
}


info=$(my_date)
echo $info
echo $1
if [$1 = ""]; then
  exit;
fi
echo "DONE"