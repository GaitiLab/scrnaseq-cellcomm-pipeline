FROM python:3.10.8-slim

RUN apt-get update && apt-get install -y procps

COPY requirements.txt /tmp/requirements.txt

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r /tmp/requirements.txt