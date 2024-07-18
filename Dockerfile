FROM nvidia/cuda:12.0.1-base-ubuntu20.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Set appropriate permissions for /tmp
RUN chmod 1777 /tmp

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    wget \
    make \
    g++ \
    gcc \
    vim \
    nano \
    libgl1-mesa-glx \
    libglib2.0-0 \
    expect \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /flask_server
COPY /flask_app_automated_docking .

RUN ulimit -s 8192
RUN wget https://boostorg.jfrog.io/artifactory/main/release/1.77.0/source/boost_1_77_0.tar.gz
RUN tar -xzf boost_1_77_0.tar.gz 
# RUN cd boost_1_77_0
# RUN ./bootstrap.sh
# RUN ./b2
# RUN ./b2 install


RUN cd /flask_server

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh \
    && bash miniconda.sh -b -p /opt/conda \
    && rm miniconda.sh

ENV PATH="/opt/conda/bin:${PATH}"
# Initialize conda for bash and create environment
RUN conda init bash \
    && conda create -y -n my_rdkit python=3.7 \
    && conda install -y -n my_rdkit -c conda-forge rdkit \
    #&& conda install -y -n my_rdkit -c conda-forge -c schrodinger pymol-bundle \
    && conda install -y -n my_rdkit -c conda-forge pymol-open-source \
    && conda install -y -n my_rdkit -c conda-forge openbabel \
    && conda run -n my_rdkit pip install meeko==0.1.dev2

# Install additional requirements
COPY requirements.txt .
RUN conda run -n my_rdkit pip install -r requirements.txt

# Copy the entrypoint script and Python script into the image
COPY entrypoint_new.sh .
# Set the entrypoint script as executable
RUN chmod +x entrypoint_new.sh


# Download the installer
RUN wget https://ccsb.scripps.edu/adfr/download/1028/ADFRsuite_Linux-x86_64_1.0_install -O /tmp/ADFRsuite_Linux-x86_64_1.0_install 

# Make the installer executable
RUN chmod a+x /tmp/ADFRsuite_Linux-x86_64_1.0_install 
# Create an expect script to handle interactive input with specific prompts
# Create an expect script to handle interactive input with specific prompts
RUN echo '#!/usr/bin/expect -f\n\
set timeout -1\n\
spawn /tmp/ADFRsuite_Linux-x86_64_1.0_install\n\
expect {\n\
    -re ".*install*" { send "Y\r"; exp_continue }\n\
    eof\n\
}\n' > /tmp/install_adfrsuite.exp

# Make the expect script executable
RUN chmod +x /tmp/install_adfrsuite.exp

# Run the expect script to install ADFRsuite
RUN /tmp/install_adfrsuite.exp

# Verify the tarball exists and extract it
RUN ls -l /flask_server/Y/ADFRsuite_x86_64Linux_1.0.tar.gz  && \
    tar zxvf /flask_server/Y/ADFRsuite_x86_64Linux_1.0.tar.gz -C /flask_server

# Verify the directory exists and run the installation script
RUN ls -l /flask_server/ADFRsuite_x86_64Linux_1.0 && \
    cd /flask_server/ADFRsuite_x86_64Linux_1.0 && \
    echo -e "y\n" | ./install.sh

# Clean up
RUN rm /tmp/ADFRsuite_Linux-x86_64_1.0_install /tmp/install_adfrsuite.exp

RUN rm -rf /flask_server/Y

#Vina
RUN wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.3/vina_1.2.3_linux_x86_64
RUN chmod a+x vina_1.2.3_linux_x86_64

# Clone QuickVina2-GPU
RUN git clone https://github.com/DeltaGroupNJUPT/QuickVina2-GPU.git


WORKDIR /flask_server



# Set the entrypoint to execute the script
ENTRYPOINT ["/flask_server/entrypoint_new.sh", "/bin/bash"]



