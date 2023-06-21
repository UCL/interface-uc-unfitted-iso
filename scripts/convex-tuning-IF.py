from netgen.geom2d import SplineGeometry, CSG2d, Circle, Rectangle
from ngsolve import *
from xfem import *
from xfem.lsetcurv import LevelSetMeshAdaptation
from math import pi

import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
#import matplotlib.pyplot as plt
plt.rc('legend',fontsize=15)
plt.rc('axes',titlesize=15)
plt.rc('axes',labelsize=15)
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
np.random.seed(123)

import sys,os
sys.path.append(os.path.realpath(''))
from unf_interf_prob import SolveIF,SolveZNoCut 
from meshes import get_geometry
from problems import diffusion_4ball,helmholtz_4ball,refsol_diffusion,refsol_Helmholtz
from GeomStudy import GeomExp
from StabTuning import StabiTuningGeom,TuningGammaGeom,StabiTuningIF
from EqualOrder import EqualOrderExp

domain_type =  "convex"
helmholtz_4ball.SetDomainType(domain_type,ref_lvl = 8)
mu = [1.0,2.0]
k = [16,2]
helmholtz_4ball.Update(mu=mu,k=k,solution=refsol_Helmholtz(mu,k))
StabiTuningIF(problem= helmholtz_4ball)





# not used

#mu = [1.0,1.0]
#mu = [10.0,2.0]
#mu = [2.0,10.0]
#mu = [2.0,20.0]
#EqualOrderExp(problem= diffusion_4ball) # for paper 
#StabiTuningGeom(problem= diffusion_4ball, order_geom= 1 )
#StabiTuningGeom(problem= diffusion_4ball, order_geom= 1 )
#TuningGammaGeom(problem= diffusion_4ball, order_geom= 1,equal_order=False )
#TuningGammaGeom(problem= diffusion_4ball, equal_order=True )
#domain_type =  "squares"
#helmholtz_4ball.SetProblemType(well_posed=True)
#helmholtz_4ball.SetDomainType(domain_type,ref_lvl = 8)
#mu = [1.0,2.0]
#k = [16,2]
#helmholtz_4ball.Update(mu=mu,k=k,solution=refsol_Helmholtz(mu,k))

