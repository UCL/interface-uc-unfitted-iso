from problems import helmholtz_4ball,refsol_Helmholtz
from EqualOrder import EqualOrderExp

# k = [1,1], mu = [2.0,2.0]
# k = [10,1], mu = [2.0,2.0] 
# k = [16,6], mu = [2.0,20.0]

# values for paper
# k = [16,1], mu = [2.0,20.0] # perfect rates


domain_type =  "squares-easy" 
helmholtz_4ball.SetDomainType(domain_type,ref_lvl = 8)
#k = [16,6]
#k = [10,2]
#k = [10,10]

ks =  [[16,1],[16,1],[16,6],[16,6]]
mus = [[2,20],[2,2],[2,20],[2,2]]

for k,mu in zip(ks,mus): 
    helmholtz_4ball.Update(mu=mu,k=k,solution=refsol_Helmholtz(mu,k))
    EqualOrderExp(problem = helmholtz_4ball,show_plots=True)  

#k = [16,6]
#mu = [2.0,20.0]
#helmholtz_4ball.Update(mu=mu,k=k,solution=refsol_Helmholtz(mu,k))
#EqualOrderExp(problem=  helmholtz_4ball,show_plots=True)  

#mu = [2.0,2.0]
#mu = [20.0,2.0]
#helmholtz_4ball.Update(mu=mu,k=k,solution=refsol_Helmholtz(mu,k))
#EqualOrderExp(problem=  helmholtz_4ball,show_plots=True)  

#mu = [20.0,2.0]
#helmholtz_4ball.Update(mu=mu,k=k,solution=refsol_Helmholtz(mu,k))
#EqualOrderExp(problem=  helmholtz_4ball,show_plots=True) 

