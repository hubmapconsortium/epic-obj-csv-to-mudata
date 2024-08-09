FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    python3 \
    python3-pip \
    python3-venv \
 && rm -rf /var/cache/apt/*

COPY docker/requirements.txt /opt/requirements.txt
RUN python3 -m venv /opt/venv \
 && . /opt/venv/bin/activate \
 && pip install -U pip setuptools wheel \
 && pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip \
 && echo '. /opt/venv/bin/activate' >> /etc/bash.bashrc

WORKDIR /opt

COPY bin/obj_feature_to_mudata.py /opt
