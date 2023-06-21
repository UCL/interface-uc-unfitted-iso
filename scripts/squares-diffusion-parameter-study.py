from problems import diffusion_4ball,refsol_diffusion
from StabTuning import TuningGammaGeom,StabiTuningIF

k = [0,0]
domain_type =  "squares" 
mu = [20.0,2.0]
diffusion_4ball.SetDomainType(domain_type,ref_lvl=8)
diffusion_4ball.Update(mu=mu,k=k,solution=refsol_diffusion(mu))

TuningGammaGeom(problem= diffusion_4ball, order_geom= 1,equal_order=False, show_plots=True ) 
TuningGammaGeom(problem= diffusion_4ball, equal_order=True, show_plots=True )
StabiTuningIF(problem= diffusion_4ball, show_plots=True )


