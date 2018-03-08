FROM centos:latest
MAINTAINER Di Wu <allwu@falcon-computing.com>

# Install necessary packages
RUN yum -y install epel-release
#RUN yum -y update
RUN yum -y install boost gflags glog java wget python python-pip
RUN pip install awscli

# Add falcon files
# do this step using a local copy of the falcon package
#RUN wget https://s3.amazonaws.com/fcs-genome-pub/release/falcon-genome-latest.tgz
#RUN tar zxf falcon-genome-latest.tgz
#RUN rm falcon-genome-latest.tgz

# add the scripts
ADD ./falcon /usr/local/falcon
ADD ./common.sh /usr/local/falcon/bin
ADD ./run.sh /usr/local/falcon/bin

ENV FALCON_HOME /usr/local/falcon
ENV PATH $FALCON_HOME/bin:$FALCON_HOME/tools/bin:$PATH
ENV LD_LIBRARY_PATH=$FALCON_HOME/lib:$LD_LIBRARY_PATH
ENV LM_LICENSE_FILE=/usr/local/falcon/license.lic
