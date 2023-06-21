FROM ngsxfem/ngsxfem:latest

USER root
RUN sed -i -re 's/([a-z]{2}\.)?archive.ubuntu.com|security.ubuntu.com/old-releases.ubuntu.com/g' /etc/apt/sources.list
RUN apt-get update && apt-get dist-upgrade -y
RUN apt-get install -y  \
  texlive-full \
  wget

WORKDIR /home/app        
RUN git clone https://github.com/UCL/interface-uc-unfitted-iso.git
        
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
        
WORKDIR /home/app/interface-uc-unfitted-iso
                
