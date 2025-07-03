# Get the base image from ubuntu
FROM ubuntu:18.04

MAINTAINER Jacob Durrant <durrantj@pitt.edu>

# Install programs through apt-get
RUN apt update
RUN apt -y install wget openbabel zip bzip2

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_23.11.0-2-Linux-x86_64.sh
RUN chmod -R a+rwx Miniconda3-py38_23.11.0-2-Linux-x86_64.sh
RUN ./Miniconda3-py38_23.11.0-2-Linux-x86_64.sh -b
RUN echo "alias python=/root/miniconda3/bin/python" >> /root/.bashrc
RUN echo "alias conda=/root/miniconda3/bin/conda" >> /root/.bashrc
RUN echo "alias pip=/root/miniconda3/bin/pip" >> /root/.bashrc

# Get python dependencies
#RUN /root/miniconda3/bin/conda install -y -c conda-forge rdkit rdkit=2020.03.1
RUN /root/miniconda3/bin/conda install -y -c conda-forge rdkit=2020.09.4
RUN /root/miniconda3/bin/conda install -y numpy=1.18.1 scipy=1.4.1
RUN /root/miniconda3/bin/pip install matplotlib==3.2.1
RUN /root/miniconda3/bin/pip install func_timeout==4.3.5

# Install mgltools source: http://mgltools.scripps.edu
RUN wget https://ccsb.scripps.edu/download/548/ -O ./mgltools_x86_64Linux2_1.5.6.tar.gz
RUN tar xvfz *.tar.gz
RUN rm *.tar.gz
# NOTE: autogrow_in_docker.py is hardcoded to /mgltools_x86_64Linux2_1.5.6
RUN cd /mgltools_x86_64Linux2_1.5.6 && ./install.sh

# Copy over autogrow source files
RUN mkdir -p /autogrow4
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
