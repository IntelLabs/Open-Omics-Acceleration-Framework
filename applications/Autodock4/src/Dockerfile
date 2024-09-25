FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget \
    vim \
    build-essential \
    libxt-dev \
    libx11-dev \
    libxmu-dev \
    libxi-dev \
    libgl-dev \
    autoconf \
    automake \
    libtool \
    csh \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /autodock4

RUN wget https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/autodocksuite-4.2.6-src.tar.gz

RUN tar -xvzf autodocksuite-4.2.6-src.tar.gz

RUN rm autodocksuite-4.2.6-src.tar.gz

WORKDIR /autodock4/src/autodock

RUN mkdir build

WORKDIR /autodock4/src/autodock/build

RUN ../configure --prefix=/autodock4_prefix

RUN make && make install

ENV PATH="/autodock4_prefix/bin:$PATH"

CMD ["bash"]
