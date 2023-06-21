FROM ngsxfem/ngsxfem:latest

USER root
RUN apt-get install -y \
  texlive-full \
  wget

WORKDIR /home/app        
RUN git clone https://github.com/UCL/interface-uc-unfitted-iso.git
        
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
        
WORKDIR /home/app/interface-uc-unfitted-iso
                
