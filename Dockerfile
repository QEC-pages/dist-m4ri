FROM ubuntu:latest

# set up timezone because the image is a new ubuntu system
ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

WORKDIR /root
RUN apt-get update \
    && apt install -y git \
    && cd \
    && git clone https://bitbucket.org/malb/m4ri.git

RUN apt install -y dh-autoreconf

RUN cd m4ri \
    && autoreconf --install \
    && ./configure --enable-openmp \
    && make \
    && make check \
    && cd .libs \
    && mv *so /usr/local/lib \
    && mv ../m4ri /usr/local/include

# clean
#RUN apt-get remove git dh-autoreconf \\
#    && rm -rf m4ri


