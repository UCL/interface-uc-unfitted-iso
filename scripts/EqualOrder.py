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
from unf_interf_prob import SolveZNoCut 

def EqualOrderExp(problem,solver="sparsecholesky",show_plots=False):
    domain_type = problem.domain_type
    problem_type = problem.problem_type  
    mu = problem.mu 
    k = problem.k
    print("Running EqualOrderExp")
   
    l2_errors_order = [] 
    grad_errors_order = [] 
    ndofs_order = [] 

    orders = [1,2,3]
    for order in orders:
        
        order_geom = order     
        l2_errors = [] 
        grad_errors = [] 
        ndofs = [] 
        n_ref_max = 8
        
        stabi_dict = None 
        if domain_type == "convex":
            stabi_dict = { }
            
            stabi_dict["gamma-CIP"] = 10**(-2)
            if order == 1:
                stabi_dict["gamma-CIP"] = 10**(-3)
            stabi_dict["gamma-GLS"] = 10**(-2)
            stabi_dict["alpha-stab"] = 1e-4
            stabi_dict["gamma-IF"] = 1e-2
            stabi_dict["gamma-data"] = 10**5
            stabi_dict["gamma-Geom"] = 1e-2

        if domain_type == "squares":
            stabi_dict = { } 

            if k[0] == k[1] == 0:

                if order == 1:
                    stabi_dict["gamma-CIP"] = 5*10**(-1)
                    stabi_dict["gamma-GLS"] = 5*10**(-1)
                elif order == 2:
                    stabi_dict["gamma-CIP"] = 5*10**(-2)
                    stabi_dict["gamma-GLS"] = 5*10**(-2)
                else:
                    stabi_dict["gamma-CIP"] = 1*10**(-3)
                    stabi_dict["gamma-GLS"] = 1*10**(-3)

                stabi_dict["alpha-stab"] = 1e-5
                stabi_dict["gamma-IF"] = 1e-3
                stabi_dict["gamma-data"] = 1e5
                stabi_dict["gamma-Geom"] = 1e-2 
                    
            else: 
                stabi_dict["alpha-stab"] = 1e-5
                stabi_dict["gamma-IF"] = 1e-2
                stabi_dict["gamma-data"] = 1e5
                stabi_dict["gamma-Geom"] = 1e-5

                if order == 3:
                    stabi_dict["gamma-CIP"] = 1*10**(-2) # 1e-3 
                    stabi_dict["gamma-GLS"] = 1*10**(-5)
                elif order == 2:
                    stabi_dict["gamma-CIP"] = 1*10**(-5) 
                    stabi_dict["gamma-GLS"] = 1*10**(-5)
                else:
                    stabi_dict["gamma-CIP"] = 1*10**(-6)
                    stabi_dict["gamma-GLS"] = 1*10**(-6)


        n_refs = n_ref_max-order
        all_refs = range(n_refs)
        for n_ref in all_refs:
            
            vtk_output = False
            if domain_type == "squares":
                if order == 3 and n_ref == all_refs[-1]:
                    vtk_output = True
            #if domain_type == "convex":
            #    if order == 3 and order_geom in [1,3] and n_ref == all_refs[-1]:
            #        vtk_output = True
            result = SolveZNoCut(problem=problem, order = order, n_refs = n_ref, order_geom=order_geom, order_dual = order, stabi_dict=stabi_dict, geom_stab_all_el = True, vtk_output = vtk_output,solver=solver)
            
            l2_errors.append(result["rel-l2-err"])
            grad_errors.append(result["rel-h1sem-err"])
            ndofs.append(result["ndof"])
 
        l2_errors_order.append(l2_errors)
        grad_errors_order.append(grad_errors)
        ndofs_order.append(ndofs)

    header_str = "ndof h rel-L2-err-B rel-H1sem-err-B"
    for idx,order in zip([0,1,2],orders):
        mesh_width = np.array(ndofs_order[idx])**(-1/2)
        name_str = problem.lset_name + "-" + problem_type + "relL2error-pq{0}".format(order)+"-mus({0},{1})".format(int(mu[0]),int(mu[1]))+"-ks({0},{1})".format(int(k[0]),int(k[1]))+".dat" 
        results = [np.array(ndofs_order[idx],dtype=float),mesh_width, np.array(l2_errors_order[idx],dtype=float), np.array(grad_errors_order[idx],dtype=float)  ]
        np.savetxt(fname ="../data/{0}".format(name_str),
                           X = np.transpose(results),
                           header = header_str,
                           comments = '') 
    if show_plots: 
        for idx,order,marker in zip( [0,1,2],orders,["o","s","+"] ): 
            plt.loglog(ndofs_order[idx], l2_errors_order[idx],label="p={0}".format(order),marker="o")
        mesh_width = np.array(ndofs_order[0])**(-1/2)
        plt.loglog(ndofs_order[0],5*mesh_width,label="$\mathcal{O}(h)$", color='gray', linestyle='dashed' )
        plt.loglog(ndofs_order[0],100*mesh_width**2,label="$\mathcal{O}(h^2)$", color='gray', linestyle='dotted')
        plt.legend()
        plt.xlabel("ndof")
        plt.show() 

