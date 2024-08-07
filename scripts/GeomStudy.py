from ngsolve import *
from xfem import *
from xfem.lsetcurv import LevelSetMeshAdaptation
from math import pi
import numpy as np
import matplotlib.pyplot as plt
plt.rc('legend',fontsize=15)
plt.rc('axes',titlesize=15)
plt.rc('axes',labelsize=15)
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
np.random.seed(123)
from unf_interf_prob import SolveZNoCut 

def GeomExp(problem, show_plots=False, skip_geom=False ):
    domain_type = problem.domain_type
    problem_type = problem.problem_type  
    mu = problem.mu 
    k = problem.k
    geom_list = [1,2,3]
    if skip_geom:
        geom_list = [1]
    print("Running GeomExp")
    for order_geom in geom_list:
        l2_errors_order = [] 
        grad_errors_order = [] 
        l2_errors_order_NoIF = [] 
        ndofs_order = [] 
        orders = [1,2,3]
        n_ref_max = 8
        
        for order in orders:
            stabi_dict = { }

            if k[1] > 2:
                stabi_dict["gamma-CIP"] = 10**(-3)
                stabi_dict["gamma-GLS"] = 10**(-3)
                stabi_dict["alpha-stab"] = 1e-4
                stabi_dict["gamma-IF"] = 1e-2
                stabi_dict["gamma-data"] = 10**5
                stabi_dict["gamma-Geom"] = 1e-4
                if order == 2:
                    #stabi_dict["gamma-CIP"] = 7.5*10**(-4)
                    #stabi_dict["gamma-GLS"] = 7.5*10**(-4)
                    stabi_dict["gamma-CIP"] = 4.5*10**(-4)
                    stabi_dict["gamma-GLS"] = 4.5*10**(-4)
                if order == 3:
                    stabi_dict["gamma-CIP"] = 4*10**(-4)
                    stabi_dict["gamma-GLS"] = 4*10**(-4)

            
                #stabi_dict["gamma-CIP"] = 10**(-3)
                #stabi_dict["gamma-GLS"] = 10**(-3)
                #stabi_dict["alpha-stab"] = 1e-4
                #stabi_dict["gamma-IF"] = 1e-2
                #stabi_dict["gamma-data"] = 10**5
                #stabi_dict["gamma-Geom"] = 1e-4
            else:
                stabi_dict["gamma-CIP"] = 10**(-2)
                if order == 1:
                    stabi_dict["gamma-CIP"] = 10**(-3)
                stabi_dict["gamma-GLS"] = 10**(-2)
                stabi_dict["alpha-stab"] = 1e-4
                stabi_dict["gamma-IF"] = 1e-2
                stabi_dict["gamma-data"] = 10**5
                stabi_dict["gamma-Geom"] = 1e-2


            l2_errors = [ ]
            grad_errors = [ ]
            ndofs = [ ]
            n_refs = n_ref_max-order
            all_refs = range(n_refs)
            for n_ref in all_refs: 
                vtk_output = False
                if skip_geom:
                    if order == 2 and n_ref == 4:
                        vtk_output = True
                    result = SolveZNoCut(problem=problem, order = order, n_refs = n_ref, order_geom=order, order_dual = order, stabi_dict=stabi_dict, geom_stab_all_el = True, vtk_output = vtk_output)
                else:
                    if order == 3 and order_geom in [1,3] and n_ref == all_refs[-1]:
                        vtk_output = True
                    result = SolveZNoCut(problem=problem, order = order, n_refs = n_ref, order_geom=order_geom, order_dual = order, stabi_dict=stabi_dict, geom_stab_all_el = True, vtk_output = vtk_output)

                l2_errors.append(result["rel-l2-err"])
                grad_errors.append(result["rel-h1sem-err"])
                ndofs.append(result["ndof"])

            l2_errors_order.append(l2_errors)
            grad_errors_order.append(grad_errors)
            ndofs_order.append(ndofs)

            mesh_width = np.array(ndofs)**(-1/2)
            if skip_geom:
                name_str = problem.lset_name + "-" + problem_type + "-p{0}".format(order)+"-q{0}".format(order)+"-mus({0},{1})".format(int(mu[0]),int(mu[1]))+"-ks({0},{1})".format(int(k[0]),int(k[1]))+".dat"
            else:
                name_str = problem.lset_name + "-" + problem_type + "-p{0}".format(order)+"-q{0}".format(order_geom)+"-mus({0},{1})".format(int(mu[0]),int(mu[1]))+"-ks({0},{1})".format(int(k[0]),int(k[1]))+".dat" 
            results = [np.array(ndofs,dtype=float),mesh_width, np.array(l2_errors,dtype=float), np.array(grad_errors,dtype=float) ]
            header_str = "ndof h rel-L2-err-B rel-H1sem-err-B"
            np.savetxt(fname ="../data/{0}".format(name_str),
                               X = np.transpose(results),
                               header = header_str,
                               comments = '')
        if show_plots:
            mesh_width = np.array(ndofs_order[0])**(-1/2)
            for idx,order,marker in zip( [0,1,2],orders,["o","s","+"] ): 
                plt.loglog(ndofs_order[idx], l2_errors_order[idx],label="p={0}".format(order),marker="o")
                #plt.loglog(ndofs_order[idx], grad_errors_order[idx],label="p={0}".format(order),marker="o",linestyle="dashed")
            plt.loglog(ndofs_order[0],5*mesh_width,label="$\mathcal{O}(h)$", color='gray', linestyle='dashed' )
            plt.loglog(ndofs_order[0],100*mesh_width**2,label="$\mathcal{O}(h^2)$", color='gray', linestyle='dotted')
            plt.legend()
            plt.xlabel("ndof")
            plt.show() 
