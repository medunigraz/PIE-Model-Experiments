FROM docker.opencarp.org/opencarp/opencarp:latest

WORKDIR /openCARP

#RUN python3 -m pip install latex

RUN rm -r examples
RUN git clone -b main https://github.com/medunigraz/PIE-Model-Experiments.git

# ADD bin ./PIE-Model-Experiments/bin
# ADD calibration ./PIE-Model-Experiments/calibration
# ADD scripts ./PIE-Model-Experiments/scripts
# ADD setups ./PIE-Model-Experiments/setups
# ADD software ./PIE-Model-Experiments/software
# ADD eval.sh ./PIE-Model-Experiments

ENV PATH="/openCARP/PIE-Model-Experiments/bin:$PATH"
#ENV PATH="/openCARP/PIE-Model-Experiments/software/ParaView/bin:$PATH"

RUN apt-get update \
    && apt-get install -y \
        nmap \
        vim




