FROM continuumio/miniconda3:23.3.1-0

RUN git clone https://github.com/iquasere/KEGGCharter.git \
&& conda install -c conda-forge mamba libarchive=3.6.2=h039dbb9_1 \
&& mamba env update --file KEGGCharter/cicd/keggcharter_env.yml --name base \
&& bash KEGGCharter/envs/ci_build.sh \
&& conda clean --all -y \
&& apt-get purge -y --auto-remove $buildDeps

CMD [ "python", "bin/keggcharter.py" ]