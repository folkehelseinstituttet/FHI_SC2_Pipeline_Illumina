FROM garcianacho/fhibaseillumina:24052022
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"

COPY CommonFiles/ /home/docker/CommonFiles/
COPY Scripts/ /home/docker/Scripts/
COPY Binaries/ /home/docker/Binaries/
COPY FHI_Gisaid/ /home/docker/FHI_Gisaid/ 
RUN mkdir /Inference \
    && mkdir /Noise \
    && mkdir -p /Models/RecombinantModel \
    && chmod +x /home/docker/Binaries/* \
    && chmod +x /home/docker/Scripts/* \
    && chmod +x /home/docker/FHI_Gisaid/* \
    && chmod -R +rwx /home/docker/CommonFiles/* \
    && chmod 777 /Inference \
    && chmod 777 /Noise \
    && chmod 777 /Models/* \
    && mv /home/docker/Binaries/* /usr/bin/
USER docker
WORKDIR /home/docker/Fastq

