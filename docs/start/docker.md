---
title: Running ZPIC from a docker container
permalink: /start/docker

layout: single

sidebar:
  nav: "docs"
---

ZPIC is also available as a Docker image that allows you to use ZPIC notebooks straightforwardly. The image is available publicly on [Docker Hub](https://hub.docker.com/repository/docker/zamb/zpic) under the name `zamb/zpic`. All of the example notebooks are included in the image.

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

## Using a different port

If you already have some server running on `localhost:8888` (e.g. another Jupyter lab session) you can change the port your new Jupyter server will be running on by changing the arguments of the `-p` option. For example, to run the server on port `9999` you would do:

```bash
docker run -p 9999:8888 -t --rm zamb/zpic
```
