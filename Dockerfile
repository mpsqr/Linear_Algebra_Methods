FROM debian:latest


RUN apt-get update && apt-get install build-essential -y && apt-get install locate -y


RUN apt-get install -y \
	libblas-dev \
	liblapacke-dev


CMD /tp/run.sh