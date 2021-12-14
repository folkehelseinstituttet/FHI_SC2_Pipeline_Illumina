FROM garcianacho/fhibaseillumina:v1
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"

COPY CommonFiles/ /home/docker/CommonFiles/
COPY Scripts/ /home/docker/Scripts/
COPY Binaries/ /home/docker/Binaries/
RUN chmod +x /home/docker/Binaries/* \
    && chmod +x /home/docker/Scripts/* \
    && chmod -R +rwx /home/docker/CommonFiles/* \
    && mv /home/docker/Binaries/* /usr/bin/
USER docker
WORKDIR /home/docker/Fastq

