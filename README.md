# interface-uc-unfitted-iso
This repository contains software and instructions to reproduce the numerical experiments in the paper
> Unique continuation for an elliptic interface problem using unfitted isoparametric finite elements
>
> * authors: Erik Burman, Janosch Preuss
> * University College London

# How to run / install
We describe two options to setup the software for running the experiments. 

* downloading a `docker image` from `Zenodo` or `Docker Hub` which contains all dependencies and tools to run the application,
* or installing everything manually on your own local machine. 

We recommend the first option as it is quick and convenient. The second option provides higher flexibility but may be more complicated. 
Please contact <j.preuss@ucl.ac.uk> if problems occur. 
 
The instructions for running the image are geared towards users who have access to a Unix like environment with a `bash` shell.
Windows users may use Linux subsystems or tools like [Git BASH](https://gitforwindows.org/) or [MobaXterm](https://mobaxterm.mobatek.net/) to 
run these commands.

## Pulling the docker image from Docker Hub 
* Please install the `docker` platform for your distribution as described [here](https://docs.docker.com/get-docker/).
* After installation the `Docker daemon` has to be started. This can either be done on boot or manually. In most Linux 
distributions the command for the latter is either `sudo systemctl start docker` or `sudo service docker start`.
* Pull the docker image using the command `docker pull janosch2888/interface-uc:v2`. 
* Run the image with `sudo docker run -it janosch2888/interface-uc:v2 bash`.
* Proceed further as described in [How to reproduce](#repro).

## Downloading the docker image from Zenodo
* For this option the first two steps are the same as above.
* The image can be downloaded [here](https://doi.org/10.5281/zenodo.8134572). 
* Assuming that `interface-uc-repro.tar` is the filename of the downloaded image, please load the image with `sudo docker load < interface-uc-repro.tar`.
* Run the image with `sudo docker run -it janosch2888/interface-uc:v2 bash`.
* Proceed further as described in [How to reproduce](#repro).

## Manual installation
We first have to install the library `ngsxfem` as described in detail [here](https://github.com/ngsxfem/ngsxfem/blob/release/INSTALLATION.md).
The results in this article have been obtained using the software at commit 210de1e407ceb074027c9b9d93c94ccb5f3526f9 (tag: v2.1.2302). This specific commit 
should be checked out when compiling from source. Alternatively, it is possible to install the software via pip at a version that is cl

The easiest option is probably installation via `pip`. Please make sure to specify the following version 

    pip install xfem==2.1.2302
 
For compiling the figures you will also need a recent `latex` distribution installed on your machine.
Now we are ready to clone the repository using 

    git clone https://github.com/UCL/interface-uc-unfitted-iso.git 

and proceed as described in [How to reproduce](#repro).

# <a name="repro"></a> How to reproduce
The `python` scripts for runnings the numerical experiments are located in the folder `scripts`.
To run an experiment we change to this folder and run the corresponding file.
After execution has finished the produced data will be available in the folder `data`.
For the purpose of comparison, the folder `data_save` contains a copy of the data which has been used for the plots in the paper.
The data in both folders should be identical.

To generate the plots as shown in the article from the data just produced we change to the folder `plots`
and compile the corresponding `latex` file.
Below we decribe the above process for each of the figures in the article in detail.
For viewing the generated pdf file, say `figure.pdf`, the figure has to be copied to the host machine.
This can be done by executing the following commands in a new terminal window (not the one in which `docker` is run):

    CONTAINER_ID=$(sudo docker ps -alq)
    sudo docker cp $CONTAINER_ID:/home/app/interface-uc-unfitted-iso/plots/figure.pdf \
    /path/on/host/machine/figure.pdf

Here, `/path/on/host/machine/` has to be adapted according to the file structure on the host machine.
The file `figure.pdf` can then be found at the designated path on the host machine and inspected with a common pdf viewer.
(The command above assumes that the reproduction image is the latest docker image to be started on the machine).
Alternatively, if a recent latex distribution is available on the host machine it is also possible to copy data and tex files to the latter and
compile the figures there.


## <a name="Fig4"></a> Figure 4
Change to directory `scripts`. Run 

    python3 squares-diffusion-contrast.py

Data files of the form "ball-4-norm-squares-iprelL2error-pq __i__ -mus(__a__,__b__)-ks(0,0).dat" will be created. Here, __i__ describes the finite element order.
The value of the diffusion parameter in subdomain Omega1 is mu1 = __a__ and in subdomain Omega2 is mu2 = __b__. The data in the files is structured in the following columns:

* ndof: Degrees of freedom in the finite element space (= CutFEM space + space for dual variable)
* h: Proportional to the width of the mesh
* rel-L2-err-B: relative L2-error in target domain B 
* rel-H1sem-err-B: relative H1-semi-error (only the gradient) in target domain B 

The raw vtk-data for the plots of the absolute error is contained in the files "ball-4-norm-squares-ip-p2-q2-mus(__a__,__b__)-ks(16,1)-lvl4.vtu", where 
mu1 = __a__ and in subdomain Omega1 is mu2 = __b__ in subdomain Omega2 (you may copy these files to your own machine and inspect them using paraview).
To generate Figure 4, switch to the folder `plots` and run  

    lualatex -pdf ball-4-norm-squares-diffusion-contrast.tex

## Figure 5
Change to directory `scripts`. Run 

    python3 squares-diffusion-parameter-study.py

This may take a while (~ 1hour). The following files will be produced:

* The file "ball-4-norm-squares-param-tuning-ip-gamma-IF-EqualOrder.dat" contains the data for the adjoint consistent method (solid line) in the left plot. 
The data for the method which is not adjoint consistent is available in "ball-4-norm-squares-param-tuning-ip-gamma-IF-EqualOrder-NoAdjCons.dat". The first column in the file is the value of the stabilization parameter, the other columns contain the relative L2-errors for the respective orders p.
* "ball-4-norm-squares-param-tuning-ip-gamma-Geom-q1.dat" contains the data for the middle plot of Figure 5, so q=1.
* "ball-4-norm-squares-param-tuning-ip-gamma-Geom-pq-equal.dat" contains the data for the rightmost plot of Figure 5

To generate Figure 5, switch to the folder `plots` and run  

    lualatex -pdf ball-4-norm-squares-diffusion-Stab.tex

## Figure 6
Change to directory `scripts`. Run 

    python3 squares-exact-data-Helmholtz.py

Data files of the form "ball-4-norm-squares-iprelL2error-pq __i__ -mus(__a__,__b__)-ks(__c__,__d__).dat". The meaning of the variables __i__,__a__ and __b__ is like 
explained in [Figure 4](#Fig4). Further,  k1 = __c__ is the wavenumber in subdomain Omega1 and k2 = __d__ in subdomain Omega2. 
Additionally, `vtu` files (called accordingly) containing the data for the plots in the inset of Figure 6 will be produced. To generate Figure 6, switch to the folder `plots` and run 

    lualatex -pdf ball-4-norm-squares-Helmholtz-contrast-rev.tex

## <a name="Fig7"></a> Figure 7
For reproduce the data, change to directory `scripts` and run

    python3 convex-exact-data.py

Afterwards, new data files of the form "ball-4-norm-convex-ip-p __i__ -q __j__ -mus(__a__,__b__)-ks(__c__,__d__).dat" will be available in the folder `data`. Here, __i__ in [1,2,3] describes the finite element order p and __j__ the order of the geometry approximation q in [1,2,3] using the same notation as in the paper. 
As above, __a__ is mu1, __b__ is mu2, __c__ is k1 and __d__ is k2.
The data in these files is ordered in the same way as described in [Figure 4](#Fig4). 

There are also `vtu` files produced (named similarly as above) which contain the raw data for the plots of the absolute error show in the insets of 
Figure 7. To generate Figure 7, switch to the folder `plots` and run 

    lualatex -pdf ball-4-norm-convex-ip-qs-pres.tex
    lualatex -pdf ball-4-norm-convex-Helmholtz-wavenumber.tex

The furst command will compile figure (A) and the second figure (C).

## Figure 8
Change to directory `scripts`. Run

    python3 convex-noise.py   

Data files of the form "ball-4-norm-convex-ip-p __i__ -q __i__ -theta __j__-deltap0__X__.dat" will be created. Here:

*  __i__ in [1,2,3] describes the finite element order p, which is chosen equal to the order of the geometry approximation q for this problem,
* __j__ in [0,1,2] gives the value of theta (as defined in Secion 5.2.2 of the paper),
* __X__ in [1,5] gives the first entry of the the variable delta_p (also defined in Section 5.2.2). So, __X__ = 1 corresponds to the solid lines 
in the plot while __X__ = 5 corresponds to the dashed lines. 

To generate Figure 8, switch to the folder `plots` and run  
    
    lualatex -pdf  ball-4-norm-convex-Helmholtz-noise.tex

## Figure 9

#### Disclaimer: 
This is a 3D experiment which is computationally intensive. The final refinement level requires the machine to have about 400Gb of RAM.
Change to directory `scripts`. Run 

    python3 3D-comp-unfitted-stab.py
    python3 3D-comp-fitted-stab.py

The data files for the unfitted case containing the relative L2-errors in the targe domain are called "3D-comp-unfitted-stab--p__j__-q__j__-mus(3,__a__)-ks(10,30).dat" where where __j__ in [1,2,3] represents the finite element order and __a__ in [3,30] the value of mu2 in subdomain Omega2. 
The data files for the fitted case containing the errors are "3D-comp-fitted-stab-ip-p__j__-q__j__-mus(3,__a__)-ks(10,30).dat". 
These are the data files for the top panel of the plot. The files for the bottom panel are also available. 
These are called "3D-comp-unfitted-stab-u-eval-reflvl2-p__j__-mus(3,__a__)-ks(10,30).dat" for the unfitted case and "3D-comp-fitted-stab-u-eval-reflvl2-p__j__-mus(3,__a__)-ks(10,30).dat" for the fitted case. The data files contain the following three columns: 

* linecord: This is the x1 coordinate which is given on the horozontal axis of the plots.
* uval: The values of the reference solution (grey in the plot).
* uhval: The values of the numerical solution.

To generate Figure 9, switch to the folder `plots` and run 

    lualatex -pdf ball-3d-fitted-unfitted-contrast.tex




