FROM continuumio/miniconda3

RUN buildDeps='build-essential zlib1g-dev' \
&& apt-get update \
&& apt-get install -y $buildDeps --no-install-recommends \
&& rm -rf /var/lib/apt/lists/* \
&& git clone https://github.com/iquasere/KEGGCharter.git \
&& conda install -c conda-forge -y mamba \
&& mamba env update --file KEGGCharter/envs/environment.yml --name base \
&& bash KEGGCharter/envs/ci_build.sh \
&& conda clean --all -y \
&& apt-get purge -y --auto-remove $buildDeps

CMD [ "python", "bin/keggcharter.py" ]