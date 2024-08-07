from problems import helmholtz_4ball,refsol_Helmholtz
from GeomStudy import GeomExp

domain_type =  "convex"
helmholtz_4ball.SetDomainType(domain_type,ref_lvl = 8)

mu = [1.0,2.0]
k = [16,2]
helmholtz_4ball.Update(mu=mu,k=k,solution=refsol_Helmholtz(mu,k))
GeomExp(problem= helmholtz_4ball,show_plots=False)

mu = [2.0,20.0]
k = [16,1]
helmholtz_4ball.Update(mu=mu,k=k,solution=refsol_Helmholtz(mu,k))
GeomExp(problem= helmholtz_4ball,show_plots=False, skip_geom=True)

mu = [2.0,20.0]
k = [16,6]
helmholtz_4ball.Update(mu=mu,k=k,solution=refsol_Helmholtz(mu,k))
GeomExp(problem= helmholtz_4ball,show_plots=False, skip_geom=True)

mu = [2.0,20.0]
k = [32,1]
helmholtz_4ball.Update(mu=mu,k=k,solution=refsol_Helmholtz(mu,k))
GeomExp(problem= helmholtz_4ball,show_plots=False, skip_geom=True)
