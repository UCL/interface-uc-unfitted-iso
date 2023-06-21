from ngsolve import *
from xfem import *
from xfem.lsetcurv import LevelSetMeshAdaptation
import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
plt.rc('legend',fontsize=15)
plt.rc('axes',titlesize=15)
plt.rc('axes',labelsize=15)
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
from unf_interf_prob import SolveZNoCut 


def TuningGammaGeom(problem,order_geom=1,equal_order=False,show_plots=False):
    
    domain_type = problem.domain_type
    problem_type = problem.problem_type  
    mu = problem.mu 
    k = problem.k
    order_geom = order_geom
    n_ref = 4
    orders = [1,2,3]
    stabi_dict = { }

    for param,val_range in zip([  "gamma-Geom"], 
                               [ [10**(-14+i) for i in range(21)] ]):

        l2_errors_order = [ ]
        l2_errors_no_adj_cons_order = [ ]
        for order in orders:
            if equal_order:
                order_geom = order
            if domain_type == "squares":
                if k[0] == k[1] == 0:
                    stabi_dict["gamma-CIP"] = 10**(-4)
                    stabi_dict["gamma-GLS"] = 10**(-4)
                    stabi_dict["alpha-stab"] = 1e-5
                    stabi_dict["gamma-IF"] = 1e0
                    stabi_dict["gamma-Geom"] = 1e-3 
                    stabi_dict["gamma-data"] = 1e5                       
                else:
                    stabi_dict["gamma-CIP"] = 5*10**(-2+2*order)
                    stabi_dict["gamma-GLS"] = 5*10**(-5+order)
                    stabi_dict["alpha-stab"] = 1e-5
                    stabi_dict["gamma-IF"] = 1e-2
                    stabi_dict["gamma-data"] = 1e5
                    stabi_dict["gamma-Geom"] = 1e-2
 
            l2_errors_param = [] 
            l2_errors_no_adj_cons = [] 
            for val in val_range:
                stabi_dict[param] = val
                result = SolveZNoCut(problem=problem, order = order, n_refs = n_ref, order_geom=order_geom, order_dual=order, stabi_dict = stabi_dict,geom_stab_all_el = True )
                l2_err = result["rel-l2-err"]
                ndof = result["ndof"]
                l2_errors_param.append( l2_err)
                if param == "gamma-IF":
                    result = SolveZNoCut(problem=problem, order = order, n_refs = n_ref, order_geom=order_geom, stabi_dict = stabi_dict,geom_stab_all_el = True, adjoint_consistent = False )       
                    l2_errors_no_adj_cons.append(result["rel-l2-err"])
            l2_errors_order.append(l2_errors_param)
            if param == "gamma-IF":
                l2_errors_no_adj_cons_order.append(l2_errors_no_adj_cons)

  
        name_str = problem.lset_name + "-param-tuning-" + problem_type + "-" + param + "-" +"q{0}".format(order_geom) + ".dat"         
        if equal_order:
            name_str = problem.lset_name + "-param-tuning-" + problem_type + "-" + param + "-" +"pq-equal" + ".dat"         

        results = [np.array( val_range ,dtype=float)] + [np.array( l2_errors_order[i-1],dtype=float) for i in orders  ]
        header_str = "param" 
        for order in orders:
            header_str += " order{0}".format(order)
        np.savetxt(fname ="../data/{0}".format(name_str),
                               X = np.transpose(results),
                               header = header_str,
                               comments = '')
         
        print("Tuning {0}".format( param ))  
        if show_plots:
            for idx,order,marker in zip( [0,1,2],orders,["o","s","+"] ):  
                plt.loglog(np.array(val_range),  np.array(l2_errors_order[idx],dtype=float) ,label="k={0}".format(order),marker="o")
            plt.legend()
            plt.xlabel("param")
            plt.title(param) 
            plt.savefig("relerr.png",transparent=True,dpi=200)
            plt.show() 
        

def StabiTuningIF(problem,order_geom=1,show_plots=False):
    
    domain_type = problem.domain_type
    problem_type = problem.problem_type  
    mu = problem.mu 
    k = problem.k
    order_geom = order_geom
    n_ref = 4
    orders = [1,2,3]
    stabi_dict = { }
    
    for param,val_range in zip([  "gamma-IF" ], 
                               [  [10**(-14+i) for i in range(28)] ]):
        l2_errors_order = [ ]
        l2_errors_no_adj_cons_order = [ ]
        for order in orders:
            order_geom = order    
            if k[0] == k[1] == 0:    
                stabi_dict["gamma-CIP"] = 10**(-4)
                stabi_dict["gamma-GLS"] = 10**(-4)
                stabi_dict["alpha-stab"] = 1e-5
                stabi_dict["gamma-IF"] = 1e-3
                stabi_dict["gamma-Geom"] = 1e-3 
                stabi_dict["gamma-data"] = 1e5
            else:
                stabi_dict["gamma-CIP"] = 5*10**(-2+2*order)
                stabi_dict["gamma-GLS"] = 5*10**(-5+order)
                stabi_dict["alpha-stab"] = 1e-5
                stabi_dict["gamma-IF"] = 1e-2
                stabi_dict["gamma-data"] = 1e5
                stabi_dict["gamma-Geom"] = 1e-2
             
            l2_errors_param = [] 
            l2_errors_no_adj_cons = [] 
            for val in val_range:
                stabi_dict[param] = val
                result = SolveZNoCut(problem=problem, order = order, n_refs = n_ref, order_geom=order_geom, order_dual=order, stabi_dict = stabi_dict,geom_stab_all_el = True )
                l2_err = result["rel-l2-err"]
                ndof = result["ndof"]
                l2_errors_param.append( l2_err)
                if param == "gamma-IF":
                    result = SolveZNoCut(problem=problem, order = order, n_refs = n_ref, order_geom=order_geom, order_dual=order, stabi_dict = stabi_dict,geom_stab_all_el = True, adjoint_consistent = False )       
                    l2_errors_no_adj_cons.append(result["rel-l2-err"])
            l2_errors_order.append(l2_errors_param)
            if param == "gamma-IF":
                l2_errors_no_adj_cons_order.append(l2_errors_no_adj_cons)
 
        name_str = problem.lset_name + "-param-tuning-" + problem_type + "-" + param + "-" +"EqualOrder" + ".dat"         
        results = [np.array( val_range ,dtype=float)] + [np.array( l2_errors_order[i-1],dtype=float) for i in orders  ]
        header_str = "param" 
        for order in orders:
            header_str += " order{0}".format(order)
        np.savetxt(fname ="../data/{0}".format(name_str),
                               X = np.transpose(results),
                               header = header_str,
                               comments = '')  
        
        if param == "gamma-IF":
            name_str = problem.lset_name + "-param-tuning-" + problem_type + "-" + param + "-" +"EqualOrder" + "-NoAdjCons.dat"         
            results = [np.array( val_range ,dtype=float)] + [np.array( l2_errors_no_adj_cons_order[i-1],dtype=float) for i in orders  ]
            header_str = "param" 
            for order in orders:
                header_str += " order{0}".format(order)
            np.savetxt(fname ="../data/{0}".format(name_str),
                                   X = np.transpose(results),
                                   header = header_str,
                                   comments = '')
            if show_plots:
                for idx,order,marker in zip( [0,1,2],orders,["o","s","+"] ):  
                    plt.loglog(np.array(val_range),  np.array(l2_errors_order[idx],dtype=float) ,label="p={0}".format(order),marker="o")
                    plt.loglog(np.array(val_range),  np.array(l2_errors_no_adj_cons_order[idx],dtype=float) ,label="p={0}".format(order),marker="o",linestyle='dashed')
                    plt.legend()
                plt.xlabel("param")
                plt.title(param) 
                plt.show() 

