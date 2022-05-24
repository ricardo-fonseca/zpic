# Docker container

ZPIC is also available as a Docker image that allows you to use ZPIC notebooks straightforwardly. The image is available publicly on [Docker Hub](https://hub.docker.com/repository/docker/zamb/zpic) under the name `zamb/zpic`.

## Launching ZPIC Jupyter notebooks from Docker

If you have Docker installed you can run this image by doing:

```bash
docker run -p 8888:8888 -t --rm zamb/zpic
```

And then launching a web browser and pointing it to `http://localhost:8888/`

## Retaining file changes

If you want to retain any file changes that happen in this session you need to mount a local directory on the docker container, e.g.

```bash
docker run -p 8888:8888 -t --rm -v $PWD:/home/jovyan/work zamb/zpic
```

All the data saved to `/home/jovyan/work` in the Docker session will be kept after the session is finished.
