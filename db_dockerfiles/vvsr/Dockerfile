FROM ubuntu:22.04

RUN  apt-get update

RUN apt-get install -y wget

RUN mkdir -p /usr/local/share/seqrepo

RUN wget --output-document=/usr/local/share/seqrepo/VV_SR_2023_05.tar https://www528.lamp.le.ac.uk/vvdata/vv_seqrepo/VV_SR_2023_05.tar

RUN tar -xvf /usr/local/share/seqrepo/VV_SR_2023_05.tar --directory /usr/local/share/seqrepo

RUN rm /usr/local/share/seqrepo/VV_SR_2023_05.tar

ENTRYPOINT []

CMD ["tail", "-f", "/dev/null"]