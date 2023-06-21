from problems import helmholtz_4ball,refsol_Helmholtz
from EqualOrder import EqualOrderExp

domain_type =  "squares" 
helmholtz_4ball.SetDomainType(domain_type,ref_lvl = 8)
k = [16,2]

mu = [2.0,20.0]
helmholtz_4ball.Update(mu=mu,k=k,solution=refsol_Helmholtz(mu,k))
EqualOrderExp(problem=  helmholtz_4ball,show_plots=True)  

mu = [2.0,2.0]
helmholtz_4ball.Update(mu=mu,k=k,solution=refsol_Helmholtz(mu,k))
EqualOrderExp(problem=  helmholtz_4ball,show_plots=True)  

mu = [20.0,2.0]
helmholtz_4ball.Update(mu=mu,k=k,solution=refsol_Helmholtz(mu,k))
EqualOrderExp(problem=  helmholtz_4ball,show_plots=True) 

