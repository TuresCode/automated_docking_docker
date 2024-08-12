# Automated Docking Docker

This project is about to create a small GPU based docker container as replacement for the bigger KAGGLE version

We try to use the tensorflow image from: https://www.tensorflow.org/install/docker

and create a Dockerfile which allows us to modify the docker container while building it.

The idea is to have an entrypoint.sh in it which runs a specific script when the docker container is started.

After creating the Dockerfile we build the docker image using this command:

docker build -t automated_docking:latest .

#start a docker container from a image
docker run -d --restart unless-stopped -it --gpus all -p 7085:7085 automated_docking

Enter the docker using:
docker exec -it docker_name /bin/bash

If the docker doesn't have access to the gpu after a while. Try to make a symlink as described here: https://stackoverflow.com/questions/72932940/failed-to-initialize-nvml-unknown-error-in-docker-after-few-hours
