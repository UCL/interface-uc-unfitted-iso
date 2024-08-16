FROM ngsxfem/ngsxfem:v2.1.2302

WORKDIR /home/app

USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

RUN apt-get update
RUN apt-get install -y software-properties-common
RUN add-apt-repository universe
RUN apt-get install -y psmisc texlive-full

WORKDIR /home/app        
RUN git clone https://github.com/UCL/interface-uc-unfitted-iso.git
        
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
        
WORKDIR /home/app/interface-uc-unfitted-iso
                
