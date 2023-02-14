# this Dockerfile was adapted from the official bakta Dockerfile
FROM alpine:3.12

RUN apk update && apk add wget tar bash \
    && wget -q -O /etc/apk/keys/sgerrand.rsa.pub https://alpine-pkgs.sgerrand.com/sgerrand.rsa.pub \
    && wget https://github.com/sgerrand/alpine-pkg-glibc/releases/download/2.32-r0/glibc-2.32-r0.apk \
    && apk add glibc-2.32-r0.apk \
    && rm glibc-2.32-r0.apk \
    && wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba \
    && touch /root/.bashrc \
    && ./bin/micromamba shell init -s bash -p /opt/conda  \
    && cp /root/.bashrc /opt/conda/bashrc

# since we're not using the Dockerfile from the bakta repo,
# we need to wget it so we have the conda environment.yml and the bakta source code
RUN wget https://github.com/oschwengers/bakta/archive/4fb88538619955e56c80eac466479403093574de.zip \
    && unzip /4fb88538619955e56c80eac466479403093574de.zip \
    && rm /4fb88538619955e56c80eac466479403093574de.zip \
    && mv /bakta-4fb88538619955e56c80eac466479403093574de /bakta

RUN echo -e '\n  - deepsig>=1.2.5' >> /bakta/environment.yml

SHELL ["bash", "-l" ,"-c"]

RUN source /opt/conda/bashrc && micromamba activate \
    && micromamba install -y -n base -f /bakta/environment.yml \
    && micromamba clean --all --yes

RUN source /opt/conda/bashrc && micromamba activate \
    && python3 -m pip install --no-cache /bakta/ \
    && rm -rf /bakta \

# the entrypoint initializes conda and then runs command line arguments
RUN echo '#!/bin/bash' > /entrypoint.sh \
    && echo 'if [[ -z $MAMBA_INITIALIZED ]]' >> /entrypoint.sh \
    && echo 'then' >> /entrypoint.sh \
    && echo 'source /opt/conda/bashrc' >> /entrypoint.sh \
    && echo 'micromamba activate' >> /entrypoint.sh \
    && echo 'export MAMBA_INITIALIZED=1' >> /entrypoint.sh \
    && echo 'fi' >> /entrypoint.sh \
    && echo '"$@"' >> /entrypoint.sh \
    && chmod +x /entrypoint.sh \

ENTRYPOINT ["/entrypoint.sh"]
