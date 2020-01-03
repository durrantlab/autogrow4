#!/bin/bash

# Copy over autogrow source code
if [ -f autogrow ]; then
    rm -r autogrow4
fi
mkdir autogrow4
cp ../RunAutogrow.py ./autogrow4/
cp -r ../autogrow/ ./autogrow4/
cp -r ../source_compounds/ ./autogrow4/
cp -r ../utility_scripts/ ./autogrow4/
cp -r ../sample_sub_scripts/ ./autogrow4/
cp -r ../tutorial/ ./autogrow4/

# Build the autogrow docker image
docker build -t autogrow4 .

# Run clean up of temporary files
rm -rf ./autogrow4/
# rm -rf ./temp_user_files/
