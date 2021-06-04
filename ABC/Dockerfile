# OS image
FROM ubuntu:18.04

# Define working directory.
WORKDIR /root

COPY . .

# Install

# set noninteractive installation
RUN export DEBIAN_FRONTEND=noninteractive
#install tzdata package
RUN apt-get update -y
RUN apt-get install -y tzdata
# set your timezone
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
RUN dpkg-reconfigure --frontend noninteractive tzdata

RUN apt-get update -y \
  && apt-get -y install python3 python-dev python3-pip\
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip
RUN pip3 install numpy pandas scipy matplotlib

RUN \
  sed -i 's/# \(.*multiverse$\)/\1/g' /etc/apt/sources.list && \
  apt-get update -qq && \
  apt-get upgrade -qq && \
  apt-get install -qq build-essential clang && \
  apt-get install -qq software-properties-common && \
  apt-get install -qq byobu curl git htop man unzip vim wget bc gawk && \
  apt-get install -qq subversion cmake && \
  apt-get install -qq libssl-dev && \
  apt-get install -qq mpich libmpich-dev && \
  apt-get install -qq awscli && \
  rm -rf /var/lib/apt/lists/*

# Set environment variables.
ENV HOME /root
# Define default command.
CMD ["bash"]

# Update gcc/g++ to v4.9
#RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test
#RUN apt-get update -qq
#RUN apt-get install -qq gcc-4.9 g++-4.9
#RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 90
#RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.9 90
#RUN update-alternatives --install /usr/bin/gcov gcov /usr/bin/gcov-4.9 90

# Install LCOV
RUN git clone https://github.com/linux-test-project/lcov.git
RUN make -C lcov install
RUN apt-get update
RUN apt-get install -y ruby-full
RUN gem install coveralls-lcov

# Update CMake to v3.5.1
RUN git clone --depth 1 -b v3.18.2 https://github.com/Kitware/CMake.git
RUN \
  cd CMake && \
  mkdir build && \
  cd build && \
  cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr && \
  make -j4 && \
  make install && \
  ldconfig && \
  cd ../.. && \
  cmake --version

# Clean working directory
RUN rm -rf openmp
RUN rm -rf CMake

# Get and install boost
RUN wget http://downloads.sourceforge.net/project/boost/boost/1.60.0/boost_1_60_0.tar.gz \
  && tar xfz boost_1_60_0.tar.gz \
  && rm boost_1_60_0.tar.gz \
  && cd boost_1_60_0 \
  && ./bootstrap.sh --with-libraries=program_options \
  && ./b2 install

# Compile and install pakman
# RUN git clone https://github.com/ThomasPak/pakman
RUN cd pakman && mkdir -p build && cd build && cmake .. && make

# Compile all programs needed to run pakman and perform ABC-SMC inference
ENV PATH=${PATH}:"pakman/build/src/"
ENV OMPI_MCA_mpi_warn_on_fork=0

RUN g++ -O3 -std=c++11 simulator.cpp population.cpp BRand.cpp -o simulator && \
  g++ -O3 -std=c++11 prior_sampler.cpp -o prior_sampler && \
  g++ -O3 -std=c++11 prior-pdf.cpp -o prior-pdf && \
  g++ -O3 -std=c++11 perturbation-pdf.cpp -o perturbation-pdf && \
  g++ -O3 -std=c++11 perturbation-pdf-ln-normal.cpp -o perturbation-pdf-ln-normal && \
  g++ -O3 -std=c++11 perturber-ln-normal.cpp -o perturber-ln-normal && \
  g++ -O3 -std=c++11 prior_sampler_x.cpp -o prior_sampler_x && \
  g++ -O3 -std=c++11 prior-pdf_x.cpp -o prior-pdf_x && \
  g++ -O3 -std=c++11 perturber_x.cpp -o perturber_x && \
  g++ -O3 -std=c++11 perturbation-pdf_x.cpp -o perturbation-pdf_x && \
  g++ -O3 -std=c++11 x_simulator_model2.cpp population_sex_diff_v2.cpp BRand.cpp -o x_simulator && \
  g++ -O3 -std=c++11 simulator_expand.cpp population.cpp BRand.cpp -o simulator_expand && \
  g++ -O3 -std=c++11 simulator_expand_smallback.cpp population.cpp BRand.cpp -o simulator_expand_smallback && \
  g++ -O3 -std=c++11 x_simulator_model2_expand.cpp population_sex_diff_v2.cpp BRand.cpp -o x_simulator_expand && \
  g++ -O3 -std=c++11 x_simulator_model2_expand_smallback.cpp population_sex_diff_v2.cpp BRand.cpp -o x_simulator_expand_smallback

# Set the env variables
ARG NAME
ENV NAME=$NAME
RUN echo "Inferring hs for gene: " $NAME

ARG CHROM
ENV CHROM=$CHROM
RUN echo "Chromosome:" $CHROM

ARG MUT
ENV MUT=$MUT
RUN echo "LoF mutation rate: " $MUT

ARG FREQ
ENV FREQ=$FREQ
RUN echo "Observed LoF frequency: " $FREQ

ARG SAMP_SIZE
ENV SAMP_SIZE=$SAMP_SIZE
RUN echo "Haploid sample size: " $SAMP_SIZE

ARG NUM
ENV NUM=$NUM
RUN echo "Observed total LoF number: " $NUM

ARG NPROC
ENV NPROC=$NPROC
RUN echo "Number of processors requested for mpi: " $NPROC

ARG NACCEPT
ENV NACCEPT=$NACCEPT
RUN echo "Number of acceptances: " $NACCEPT

RUN chmod +x run_abc.sh
RUN mkdir -p outfiles

#ADD fetch_and_run.sh /usr/local/bin/fetch_and_run.sh

RUN chmod +x fetch_and_run.sh

ENTRYPOINT ["./fetch_and_run.sh"]

#CMD ["sh", "-c", "./run_abc.sh $NAME $CHROM $MUT $FREQ $SAMP_SIZE $NUM $NPROC $NACCEPT"]
