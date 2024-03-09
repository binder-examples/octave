FROM quay.io/jupyterhub/repo2docker:2023.06.0-75.g3221560@sha256:3f0a789f7f779ebd4d57fa62e8fc669b4761cf5b51ece0cabe9494257f5ac6e9

RUN apk add --no-cache curl build-base python3 python3-dev py3-pip

RUN python3 -m pip install --upgrade wheel setuptools

# https://stackoverflow.com/a/41651363/1695486
RUN apk add --no-cache curl curl-dev
COPY create_docker_image.sh /create_docker_image.sh
COPY binder_cache.py /binder_cache.py
COPY trigger_binder.sh /trigger_binder.sh

ENTRYPOINT ["/bin/bash", "/create_docker_image.sh"]
