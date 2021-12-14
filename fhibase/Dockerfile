FROM debian:bullseye
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"
ENV PATH="/home/docker/miniconda3/envs/pangolin/bin:/home/docker/Scripts:/home/docker/miniconda3/bin:${PATH}"
ARG PATH="/home/docker/miniconda3/envs/pangolin/bin:/home/docker/Scripts:/home/docker/miniconda3/bin:${PATH}"

#Based on r-base 4.05's Dockerfile
RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& addgroup docker staff

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
		fonts-texgyre \
        liblzma-dev \
        libcurl4-openssl-dev \
        libxml2-dev \   
        libgc-dev \
        libssl-dev \
        gcc-9 \
        gcc-9-base \
        libpoppler-cpp-dev \
        nano \
        procps \
        git-all \
        bc \
        rsync \
        nodejs \
        npm \
    	&& rm -rf /var/lib/apt/lists/*

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

ENV R_BASE_VERSION 4.0.5

## During the freeze, new (source) packages are in experimental and we place the binaries in our PPA
RUN echo "deb http://deb.debian.org/debian experimental main" > /etc/apt/sources.list.d/experimental.list \
    && echo "deb [trusted=yes] https://eddelbuettel.github.io/ppaR400 ./" > /etc/apt/sources.list.d/edd-r4.list

## Now install R and littler, and create a link for littler in /usr/local/bin
RUN apt-get update \
        && apt-get install -y --no-install-recommends \
                gcc-9-base \
        libopenblas0-pthread \
		littler \
        r-cran-littler \
		r-base \
		r-base-dev \
        r-base-core\
		r-recommended \
	&& ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installBioc.r /usr/local/bin/installBioc.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installDeps.r /usr/local/bin/installDeps.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
	&& install.r docopt \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& rm -rf /var/lib/apt/lists/*

RUN npm install --global @nextstrain/nextclade@0.14.4 \
	&& userdel docker && useradd -ms /bin/bash docker

ARG USER=docker
ARG GROUP=docker
ARG UID
ARG GID

ENV USER=$USER
ENV GROUP=$GROUP
ENV UID=$UID
ENV GID=$GID
ENV HOME="/home/${USER}"
