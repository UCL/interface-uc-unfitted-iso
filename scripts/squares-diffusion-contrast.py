from problems import diffusion_4ball,refsol_diffusion
from EqualOrder import EqualOrderExp

k = [0,0]
domain_type =  "squares" 

mu = [2.0,2.0]
diffusion_4ball.SetDomainType(domain_type,ref_lvl=8)
diffusion_4ball.Update(mu=mu,k=k,solution=refsol_diffusion(mu))
EqualOrderExp(problem= diffusion_4ball)  

mu = [2.0,20.0]
diffusion_4ball.SetDomainType(domain_type,ref_lvl=8)
diffusion_4ball.Update(mu=mu,k=k,solution=refsol_diffusion(mu))
EqualOrderExp(problem= diffusion_4ball) # for paper 

mu = [20.0,2.0]
diffusion_4ball.SetDomainType(domain_type,ref_lvl=8)
diffusion_4ball.Update(mu=mu,k=k,solution=refsol_diffusion(mu))
EqualOrderExp(problem= diffusion_4ball) # for paper 



