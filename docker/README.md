This directory contains the scripts to run AutoGrow4 within a docker container.
This is useful when using OS that is not compatible with AutoGrow dependencies or
are not supported by our current multithreading approach, such as Windows.

Prior to running these scripts please install docker software. Please also
always run these scripts with sudo (Linux/MacOS) or administrative
priveledges (Window).

The first time running AutoGrow4 with docker (or if docker images have been purged),
will take a few minutes longer to install dependencies within the docker environment.

Depending on the AutoGrow4 settings, processor speed/count... AutoGrow4 may complete
within minutes or may take as long as multiple days. For this reason please make sure
to run in settings that can handle that will be able to run complete without being
disrupted or disruptive. Using nohup may be a useful wrapper for longer runs or
when running jobs remotes (ie running a job over ssh).

# Run instructions
To run AutoGrow4 in a docker, please run the `autogrow_in_docker.py` script:
    Example on Linux/MacOS:
        #  cd to this directory in a bash terminal
        1) cd autogrow4/docker/
        # Run autogrow_in_docker.py with sudo and supply a json file using the
        # normal pathing of your system.
        # Please note that the docker downloads its own copy of obabel and MGLTools
        # so you do not need to provide those paths.
        2) `sudo python autogrow_in_docker.py -j ./examples/sample_autogrow_docker_json.json`

        # Results will be output to the directory specified by the root_output_folder variable

    Example on Windows OS:
        1) open a docker enabled and bash enabled terminal with administrative priveledges
        #  cd to this directory in a bash terminal
        3) cd autogrow4/docker/
        4)  `python autogrow_in_docker.py -j ./examples/sample_autogrow_docker_json.json`

        # Results will be output to the directory specified by the root_output_folder variable

Files
=====

For Use in the Host System
--------------------------

* `autogrow_in_docker.py`: Run AutoGrow from within docker. Launches docker
  image. Accepts the exact same parameters as AutoGrow4, with the following
  exceptions:
    1) User variables must be supplied in JSON format.
        - Please see documentation within the tutorial manual and an example can be found:
          -  ./examples/sample_autogrow_docker_json.json

    Required variables within the JSON file:
    - `-root_output_folder`: folder path on host system that results will be copied to.
    - `-source_compound_file`: Path on host system to the tab-delineate .smi file that will seed generation 1.
    - `-filename_of_receptor`: Path on host system of the receptor to be tested.
    - `-center_x`, `-center_y`, `-center_z`: x,y,z coordinates of center of pocket to be tested.
    - `-size_x`, `-size_y`, `-size_z`: dimensions of the pocket in x,y,z coordinates.
    Variable that will be ignored:
    - `-openbabel_bin_directory` should not be specified.
    - `-mgltools_directory` should not be specified.

* `examples/example.bash`: An example of how to run `autogrow_in_docker.py`.
* `examples/sample_autogrow_docker_json.json`: A sample JSON file to supply `autogrow_in_docker.py`.

For Use in Docker
-----------------

* `run_autogrow_in_container.bash`: The docker image's ENTRYPOINT runs this script.
* `run_autogrow_in_container_windows.bash`: The windows version of docker image's ENTRYPOINT runs this script. It is automatically switched by autogrow_in_docker.py
* `Dockerfile`: Docker instructions re. how to build the image.