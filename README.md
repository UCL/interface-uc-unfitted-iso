# interface-uc-unfitted-iso
This repository contains software and instructions to reproduce the numerical experiments in the paper
> Unique continuation for an elliptic interface problem using unfitted isoparametric finite elements
>
> * authors: Erik Burman, Janosch Preuss
> * University College London

# How to run / install

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

## <a name="Fig3"></a> Fig. 3
Change to directory `scripts`. Run

    python3 convex-exact-data.py

Afterwards, new data files of the form "ball-4-norm-convex-ip-p__i__-q__j__-mus(1,2)-ks(16,2).dat" will be available in the folder `data`. Here, __i__ in [1,2,3] describes the finite element order p and __j__ the order of the geometry approximation q in [1,2,3] using the same notation as in the paper. 
The data in the files is structured in the following columns:

* ndof: Degrees of freedom in the finite element space (= CutFEM space + space for dual variable)
* h: Proportional to the width of the mesh
* rel-L2-err-B: relative L2-error in target domain B 
* rel-H1sem-err-B: relative H1-semi-error (only the gradient) in target domain B 

There are also two `vtu` files produced (named similarly as above) which contain the raw data for the plots of the absolute error show in the insets of 
Fig 3 (you may copy these files to your local machine and inspect using paraview).
To generate Fig 3, switch to the folder `plots` and run 

    lualatex -pdf ball-4-norm-convex-ip-qs-pres.tex


## Fig. 4
Change to directory `scripts`. Run

    python3 convex-noise.py   

Data files of the form "ball-4-norm-convex-ip-p__i__-q__i__-theta__j__-deltap0__X__.dat". Here:

*  __i__ in [1,2,3] describes the finite element order p, which is chosen equal to the order of the geometry approximation q for this problem,
* __j__ in [0,1,2] gives the value of theta (as defined in Secion 5.1.2 of the paper),
* __X__ in [1,5] gives the first entry of the the variable delta_p (also defined in Section 5.1.2). So, __X__ = 1 corresponds to the solid lines 
in the plot while __X__ = 5 corresponds to the dashed lines. 

To generate Fig 4, switch to the folder `plots` and run  
    
    lualatex -pdf  ball-4-norm-convex-Helmholtz-noise.tex

## <a name="Fig6"></a> Fig. 6
Change to directory `scripts`. Run 

    python3 squares-diffusion-contrast.py

Data files of the form "ball-4-norm-squares-iprelL2error-pq__i__-mus(__a__,__b__)-ks(0,0).dat" will be created. Here, __i__ describes the finite element order.
The value of the diffusion parameter in subdomain Omega1 is mu1 = __a__ and in subdomain Omega2 is mu2 = __b__. The data in these files is ordered in the 
same way as described in [Fig.3](#Fig3)
To generate Fig 6, switch to the folder `plots` and run  

    lualatex -pdf ball-4-norm-squares-diffusion-contrast.tex

## Fig. 7
Change to directory `scripts`. Run 

    python3 squares-diffusion-parameter-study.py

This may take a while (~ 1hour). The following files will be produced:

* The file "ball-4-norm-squares-param-tuning-ip-gamma-IF-EqualOrder.dat" contains the data for the adjoint consistent method (solid line) in the left plot. 
The data for the method which is not adjoint consistent is available in "ball-4-norm-squares-param-tuning-ip-gamma-IF-EqualOrder-NoAdjCons.dat". The first column in the file is the value of the stabilization parameter, the other columns contain the relative L2-errors for the respective orders p.
* "ball-4-norm-squares-param-tuning-ip-gamma-Geom-q1.dat" contains the data for the middle plot of Fig. 7, so q=1.
* "ball-4-norm-squares-param-tuning-ip-gamma-Geom-pq-equal.dat" contains the data for the rightmost plot of Fig. 7

To generate Fig 7, switch to the folder `plots` and run  

    lualatex -pdf ball-4-norm-squares-diffusion-Stab.tex


## Fig. 8
Change to directory `scripts`. Run 

    python3 squares-exact-data-Helmholtz.py

Data files of the form "ball-4-norm-squares-iprelL2error-pq__i__-mus(__a__,__b__)-ks(16,2).dat". The meaning of the variables __i__,__a__ and __b__ is like 
explained in [Fig.6](#Fig6). To generate Fig 8, switch to the folder `plots` and run  

    lualatex -pdf ball-4-norm-squares-diffusion-contrast.tex

## Fig. 9

#### Disclaimer: 
This is a 3D experiment which is computationally intensive. Please do not attempt to run this experiment 
on machines that have less than 64GB of RAM.

Change to directory `scripts`. Run 

    python3 3DStudy.py

Data files of the form "ball-4-norm-convex-3D-ip-p__j__-q__j__-mus(1,2)-ks(0,0).dat" will be produced where __j__ in [1,2,3] represents the finite
element order. To generate Fig. 9, switch to the folder `plots` and run 

    lualatex -pdf ball-4-norm-convex-3D-diffusion.tex 



