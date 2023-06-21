from problems import helmholtz_4ball,refsol_Helmholtz
from NoiseExp import NoiseExp

domain_type =  "convex"
helmholtz_4ball.SetDomainType(domain_type,ref_lvl = 8)
mu = [1.0,2.0]
k = [16,2]
helmholtz_4ball.Update(mu=mu,k=k,solution=refsol_Helmholtz(mu,k))
NoiseExp(problem= helmholtz_4ball,delta_p=[1,8,28],show_plots=True)
NoiseExp(problem= helmholtz_4ball,delta_p=[5,24,80])

