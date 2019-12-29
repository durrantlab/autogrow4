Files
=====

For Use in the Host System
--------------------------

* `build.bash`: Builds the `autogrow` docker image.
* `example.bash`: An example of how to use `autogrow_in_docker.py`.
* `autogrow_in_docker.py`: Run AutoGrow from within docker. Launches docker
  image. Accepts the exact same parameters as AutoGrow, with the following
  exceptions:
  - `-vina_executable` should not be specified.
  - `-openbabel_bin_directory` should not be specified.
  - `-mgltools_directory` should not be specified.
  - `-output_dir` should not be specified. The output will be saved to
    `autogrow_output.zip` in the directory from which the
    `autogrow_in_docker.py` script is run.
  - `-directory_of_fragments` should be specified and should indicate a
    directory containing molecular fragments. Alternatively, if the user
    specifies `MW_150`, `MW_200`, or `MW_250`, the docker image will use the
    corresponding default AutoGrow fragment library (recommended).

For Use in Docker
-----------------

* `run_autogrow_in_container.bash`: The docker image's ENTRYPOINT runs this script.
* `Dockerfile`: Docker instructions re. how to build the image.