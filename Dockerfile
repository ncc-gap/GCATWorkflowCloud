FROM python:3.8.3-buster
MAINTAINER aokad <aokada@ncc.go.jp>

WORKDIR /tools
RUN apt-get -y update && \
    apt-get install -y wget unzip git

RUN wget https://github.com/aokad/ecsub/archive/v0.0.21.zip && \
    unzip v0.0.21.zip && \
    rm v0.0.21.zip && \
    cd ecsub-0.0.21; python setup.py build install

RUN git clone https://github.com/ncc-ccat-gap/GCATWorkflowCloud.git && \
    cd GCATWorkflowCloud; python setup.py build install

CMD ["/bin/bash"]
