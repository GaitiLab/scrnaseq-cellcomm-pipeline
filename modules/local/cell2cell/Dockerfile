FROM python:3.10.8-slim

COPY requirements.txt /tmp/requirements.txt

RUN apt-get update && \
    # procps as needed for Nextflow
    apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/* && \
    pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r /tmp/requirements.txt