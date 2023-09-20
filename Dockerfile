FROM continuumio/miniconda3:23.3.1-0

RUN git clone https://github.com/iquasere/KEGGCharter.git \
&& conda env update --file KEGGCharter/envs/mamba_env.yml --name base \
&& mamba env update --file KEGGCharter/envs/keggcharter_env.yml --name base \
&& bash KEGGCharter/envs/ci_build.sh \
&& conda clean --all -y \
&& apt-get purge -y --auto-remove $buildDeps

CMD [ "python", "bin/keggcharter.py" ]