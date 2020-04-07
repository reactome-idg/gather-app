FROM ubuntu:bionic
ENV DEBIAN_FRONTEND=noninteractive 

RUN apt-get update && apt-get install -y gnupg2 software-properties-common && \ 
    apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

RUN apt-get update && apt-get install -y --no-install-recommends build-essential \
    wget \
    vim \
    r-base-core \
    r-base-dev \
    python3 \
    python3-pip \
    python3-setuptools \
    python3-dev \
    pandoc \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --upgrade pip
RUN pip install pandas numpy ipython argparse reactome2py scipy statsmodels --upgrade
RUN ln -s /usr/bin/python3 /usr/bin/python && \ 
    ln -s /usr/bin/pip3 /usr/bin/pip

COPY rdependency/installPackages.R /opt/installPackages.R
RUN Rscript /opt/installPackages.R
RUN rm /opt/installPackages.R

ARG MCL_VERSION='14-137'

RUN mkdir -p /tmp/mcl && \
	cd /tmp/mcl && \
	wget https://micans.org/mcl/src/mcl-$MCL_VERSION.tar.gz && \
	tar zxf mcl-$MCL_VERSION.tar.gz && \
	cd mcl-$MCL_VERSION/ && \
	env CPPFLAGS=-I/usr/local/include LDFLAGS=-s ./configure && \
	make && \ 
	make install && \
	rm -rf /tmp/mcl  

WORKDIR /app

ENTRYPOINT ["/bin/bash"]
