FROM ubuntu:22.04

RUN  apt-get update

RUN apt-get install -y wget

RUN mkdir -p /usr/local/share/seqdata

RUN wget --output-document=/usr/local/share/seqdata/VV_SR_2024_01.tar https://www528.lamp.le.ac.uk/vvdata/vv_seqrepo/VV_SR_2024_01.tar

RUN tar -xvf /usr/local/share/seqdata/VV_SR_2024_01.tar --directory /usr/local/share/seqdata

RUN rm /usr/local/share/seqdata/VV_SR_2024_01.tar

ENTRYPOINT []

CMD ["tail", "-f", "/dev/null"]