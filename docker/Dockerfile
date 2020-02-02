# Get the base image from ubuntu
FROM ubuntu:18.04

MAINTAINER Jacob Durrant <durrantj@pitt.edu>

# Install programs through apt-get
RUN apt-get -yqq update
RUN apt-get -fyqq install
RUN apt-get -yqq install wget
RUN apt-get -yqq update
RUN apt-get -yqq install openbabel
RUN apt-get -yqq install zip
# RUN apt-get -yqq install git

# Install python
RUN apt-get -yqq install bzip2
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN chmod -R a+rwx Miniconda3-latest-Linux-x86_64.sh
RUN ./Miniconda3-latest-Linux-x86_64.sh -b
RUN rm ./Miniconda3-latest-Linux-x86_64.sh
RUN echo "alias python=/root/miniconda3/bin/python" >> /root/.bashrc


# Get python dependencies
RUN /root/miniconda3/bin/conda install -y -c conda-forge rdkit
RUN /root/miniconda3/bin/conda install -y numpy
RUN /root/miniconda3/bin/conda install -y scipy
RUN /root/miniconda3/bin/pip install matplotlib
RUN /root/miniconda3/bin/pip install func-timeout


# Install mgltools
RUN wget http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_x86_64Linux2_1.5.6.tar.gz
RUN tar xvfz *.tar.gz
RUN rm *.tar.gz
RUN cd /mgltool*/ && \
    ./install.sh
# RUN echo "alias python=/mgltools*/bin/pythonsh" >> /root/.bashrc

# Copy over autogrow files
RUN mkdir /autogrow4/
ADD ./autogrow4 /autogrow4

# Copy over User files for AutoGrow Run
RUN mkdir /UserFiles
ADD ./temp_user_files /UserFiles

# Copy over the autogrow run script
ADD ./run_autogrow_in_container.bash /autogrow/run_autogrow_in_container.bash
RUN chmod -R a+rwx /autogrow/run_autogrow_in_container.bash
ADD ./run_autogrow_in_container_windows.bash /autogrow/run_autogrow_in_container_windows.bash
RUN chmod -R a+rwx /autogrow/run_autogrow_in_container_windows.bash

# Make directories for mounted host file systems
RUN mkdir /Outputfolder


RUN ls autogrow4
ENTRYPOINT ["bash", "/autogrow/run_autogrow_in_container.bash"]

# Windows version. Automatically turned on by autogrow_in_docker.py
# ENTRYPOINT ["bash", "/autogrow/run_autogrow_in_container_windows.bash"]


# ENTRYPOINT ["bash"]