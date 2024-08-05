from ngsolve import *
from xfem import *
from xfem.lsetcurv import LevelSetMeshAdaptation
import numpy as np
import matplotlib.pyplot as plt
plt.rc('legend',fontsize=15)
plt.rc('axes',titlesize=15)
plt.rc('axes',labelsize=15)
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
np.random.seed(123)

from unf_interf_prob import SolveZNoCut,SolveFitted 
from meshes import get_geometry
from problems import interface_problem, levelset_2ball, refsol_Helmholtz_2ball

def SolvedFittedConcentric(problem,show_plots=False,r_eval=[]):
    domain_type = problem.domain_type
    problem_type = problem.problem_type  
    mu = problem.mu 
    k = problem.k
    print("Running ResolvedGeom")
    l2_errors_order = [] 
    h1s_errors_order = [] 
    l2_errors_order_NoIF = [] 
    ndofs_order = []  
    #orders = [1,2,3]
    orders = [1]
    #n_ref_max = 5
    n_ref_max = 4    
    mat2idx  = {'rest': 1,
                'B_outer': 1,
                'B_inner': 0,
                'omega': 0
                }
    for order in orders:
        order_geom = order
        stabi_dict = { }
        stabi_dict["gamma-CIP"] = 1e-3
        stabi_dict["gamma-GLS"] = 1e-3
        if order == 3:
            stabi_dict["gamma-CIP"] = 1e-3
            stabi_dict["gamma-GLS"] = 1e-3
        stabi_dict["alpha-stab"] = 1e-3
        stabi_dict["gamma-IF"] = 1e-4
        stabi_dict["gamma-data"] = 1e5
        stabi_dict["gamma-Geom"] = 1e-2
        
        l2_errors = [ ]
        h1s_errors = [] 
        ndofs = [ ]
        n_refs = n_ref_max-order
        vtk_output = False
        for n_ref in range(n_refs):
            #if order == 3 and n_ref ==  (n_refs-1):
            #    vtk_output = True
            #if order == 1 and n_ref == (n_refs-1):
            #   vtk_output = True
            if order == 1 and n_ref == 0:
                vtk_output = True
            result = SolveFitted(problem=problem, order = order, n_refs = n_ref, order_geom=order_geom, order_dual = order, stabi_dict=stabi_dict, geom_stab_all_el = False,vtk_output = vtk_output,mat2idx = mat2idx)
            l2_err = result["rel-l2-err"]
            ndof = result["ndof"]
            l2_errors.append( l2_err)
            ndofs.append(ndof)
            h1s_errors.append(result["rel-h1sem-err"])
            input("")
            #plt.plot(r_eval, result["u-vals"], label= "u(r)")
            #plt.plot(r_eval, result["uh-vals"], label= "uh(r)")
            #plt.show() 

            results_eval = [np.array(r_eval,dtype=float), np.array(result["u-vals"],dtype=float), np.array(result["uh-vals"],dtype=float)] 
            header_str_eval = "linecord uval uhval"
            name_str_eval = "3D-comp-fitted-stab-u-eval" + "-" + "reflvl{0}".format(n_ref)+"-p{0}".format(order)+"-mus({0},{1})".format(int(mu[0]),int(mu[1]))+"-ks({0},{1})".format(int(k[0]),int(k[1]))+".dat" 
            np.savetxt(fname ="../data/{0}".format(name_str_eval),
                           X = np.transpose(results_eval),
                           header = header_str_eval,
                           comments = '')
        

        l2_errors_order.append(l2_errors)
        h1s_errors_order.append(h1s_errors) 
        ndofs_order.append(ndofs)
        print("ndofs = ", ndofs)
        print(" l2_errors = ", l2_errors )
        eoc = [log(l2_errors[i-1]/l2_errors[i])/log(2) for i in range(1,len(l2_errors))] 
        print(" eoc = ", eoc) 

        mesh_width = np.array(ndofs)**(-1/3)
        name_str = "3D-comp-fitted-stab" + "-" + problem_type + "-p{0}".format(order)+"-q{0}".format(order_geom)+"-mus({0},{1})".format(int(mu[0]),int(mu[1]))+"-ks({0},{1})".format(int(k[0]),int(k[1]))+".dat" 
        results = [np.array(ndofs,dtype=float),mesh_width, np.array(l2_errors,dtype=float), np.array(h1s_errors,dtype=float)]
        header_str = "ndof h rel-L2-err-B rel-H1s-err-B"
        np.savetxt(fname ="../data/{0}".format(name_str),
                           X = np.transpose(results),
                           header = header_str,
                           comments = '')
    if show_plots:
        mesh_width = np.array(ndofs_order[0])**(-1/problem.meshes[0].dim)
        #for idx,order,marker in zip( [0,1,2],orders,["o","s","+"] ): 
        #    plt.loglog(ndofs_order[idx], l2_errors_order[idx],label="k={0}".format(order),marker="o")
        plt.loglog(ndofs_order[0], l2_errors_order[0],label="k={0}".format(1),marker="o")
        plt.loglog(ndofs_order[0],5*mesh_width,label="$\mathcal{O}(h)$", color='gray', linestyle='dashed' )
        plt.loglog(ndofs_order[0],100*mesh_width**2,label="$\mathcal{O}(h^2)$", color='gray', linestyle='dotted')
        plt.legend()
        plt.xlabel("ndof")
        plt.show() 


domain_type = "concentric-3D-fitted"

mu = [3.0,30.0]
k = [10,30]

helmholtz_3D_ball = interface_problem(lset = levelset_2ball,
                                    solution = refsol_Helmholtz_2ball(mu,k),
                                    mu = mu,
                                    k = k,
                                    dim=3
                                      )
helmholtz_3D_ball.SetProblemType(well_posed=False)
helmholtz_3D_ball.SetDomainType(domain_type,ref_lvl = 3)

N_total = 1000
r_eval = [i*1.5/N_total for i in range(N_total-1)]
helmholtz_3D_ball.eval_pts= [ (rr, 0.0, 0.0 ) for rr in r_eval ]

#mu = [1.0,10.0]
#k = [1.0,1.0]
helmholtz_3D_ball.Update(mu=mu,k=k,solution= refsol_Helmholtz_2ball(mu=mu,k=k))
SolvedFittedConcentric(problem=helmholtz_3D_ball ,show_plots=True, r_eval=r_eval)
