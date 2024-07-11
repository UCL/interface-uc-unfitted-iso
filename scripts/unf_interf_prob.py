from ngsolve import *
from xfem import *
from xfem.lsetcurv import LevelSetMeshAdaptation
import numpy as np
np.random.seed(123)

def SolveZNoCut(problem, order = 1, n_refs = 0, order_geom = 1, theta_perturb = None, order_dual = 1, stabi_dict = None, geom_stab_all_el = True, vtk_output = False,adjoint_consistent=True,
                delta_p=None,solver="sparsecholesky" ): 

    # retrieving stabilization parameters
    if stabi_dict != None:
        gamma_CIP_order = stabi_dict["gamma-CIP"]
        gamma_GLS_order = stabi_dict["gamma-GLS"]
        alpha_stab = stabi_dict["alpha-stab"]
        gamma_IF = stabi_dict["gamma-IF"]
        gamma_data = stabi_dict["gamma-data"]
        gamma_geom = stabi_dict["gamma-Geom"] 
   
    # retrieve problem parameters
    mu = problem.mu
    mubar = (mu[0]+mu[1])/2 
    k = problem.k
    kappa = ( mu[1]/sum(mu) , mu[0]/sum(mu) )
    solution = problem.solution
    sol_gradient = problem.gradient
    coef_f = problem.coef_f
    mesh = problem.meshes[n_refs]

    #if problem.domain_type == "squares":
    #    domain_values = {'only_B': 1.0, 'omega':1.0,'full':0.0}
    #elif problem.domain_type in ["convex", "data-half", "convex-3D", "non-convex-3D"]:
    #    domain_values = {'only_B': 1.0, 'omega':1.0,'rest':0.0}
    #elif problem.domain_type == "data-all-around": 
    #    domain_values = {'only_B': 0.0, 'omega':1.0}
    #else:
    #    domain_values = {'only_B': 1.0, 'omega':1.0}
    #values_list = [domain_values[mat]
    #                   for mat in mesh.GetMaterials()]
    #cf = CoefficientFunction(values_list)

    # active mesh trafo for higher order
    if order_geom > 1:
        lsetadap = LevelSetMeshAdaptation(mesh, order=order_geom, levelset=problem.lset)
    else:
        print("Using piecewise linear reference geometry")
        lsetadap = NoDeformation(mesh, problem.lset)
    lsetp1 = lsetadap.lset_p1


    # background FE space for primal variable
    if problem.well_posed:
        Vh = H1(mesh, order=order, dirichlet="bc_Omega",dgjumps=True)
    else:
        Vh = H1(mesh, order=order, dirichlet=[],dgjumps=True)
    
    # continuous space for dual variable
    Vh0 = H1(mesh, order=order_dual, dirichlet="bc_Omega",dgjumps=False)

    # Gathering information on cut elements:
    ci = CutInfo(mesh, lsetp1)

    # Get Cut info and setup CutFEM space for primal vaiable
    hasneg = ci.GetElementsOfType(HASNEG)
    haspos = ci.GetElementsOfType(HASPOS)
    hasif = ci.GetElementsOfType(IF)
    Vh_Gamma = Compress(Vh, GetDofsOfElements(Vh, hasneg)) \
              * Compress(Vh, GetDofsOfElements(Vh, haspos)) \
              * Vh0

    # for storing solution
    gfu = GridFunction(Vh_Gamma)

    # test and trial functions
    u1,u2,z0 =  Vh_Gamma.TrialFunction()
    v1,v2,w0 =  Vh_Gamma.TestFunction()
    u = [u1,u2]
    z = [z0,z0]
    v = [v1,v2]
    w = [w0,w0]
    gradu, gradz, gradv, gradw = [[grad(fun[i]) for i in [0, 1]] for fun in [u, z, v, w]]
    gfuh = gfu.components

    # Coefficients / parameters:
    n = 1.0 / grad(lsetp1).Norm() * grad(lsetp1)
    h = specialcf.mesh_size
    nF = specialcf.normal(mesh.dim)
    # Nitsche stabilization parameter:
    lambda_nitsche  = 20
    stab = lambda_nitsche * (mu[1] + mu[0]) / h

    def P(fun):
        return fun - (fun * n) * n
    
    # jump operators across IF    
    average_flux_z = sum([- kappa[i] * mu[i] * gradz[i] * n for i in [0, 1]])
    average_flux_w = sum([- kappa[i] * mu[i] * gradw[i] * n for i in [0, 1]])
    jump_flux_u =  (mu[0] * gradu[0] - mu[1] * gradu[1]) * n 
    jump_flux_v =  (mu[0] * gradv[0] - mu[1] * gradv[1]) * n 
    jump_tangential_u =  P(gradu[0]) - P(gradu[1])
    jump_tangential_v =  P(gradv[0]) - P(gradv[1])

    # differential operator for GLS
    def calL(fun):
        hesse = [fun[i].Operator("hesse") for i in [0,1]]
        if  mesh.dim == 2:
            return (-mu[0]*hesse[0][0,0]-mu[0]*hesse[0][1,1] - k[0]**2*fun[0] ,-mu[1]*hesse[1][0,0]-mu[1]*hesse[1][1,1] - k[1]**2*fun[1]  )
        else:
            return (-mu[0]*hesse[0][0,0]-mu[0]*hesse[0][1,1]-mu[0]*hesse[0][2,2] - k[0]**2*fun[0] ,-mu[1]*hesse[1][0,0]-mu[1]*hesse[1][1,1] -mu[1]*hesse[1][2,2] - k[1]**2*fun[1]  )

    # measures for cut integrals
    dC = tuple([dCut(lsetp1, dt, deformation=lsetadap.deform,
                     definedonelements=ci.GetElementsOfType(HAS(dt)))
                for dt in [NEG, POS]])
    ds = dCut(lsetp1, IF, deformation=lsetadap.deform)

    ba_hasneg_facets = GetFacetsWithNeighborTypes(mesh, a=hasneg, b=hasneg)
    ba_haspos_facets = GetFacetsWithNeighborTypes(mesh, a=haspos, b=haspos)
    ba_facets = { NEG: ba_hasneg_facets, 
                  POS: ba_haspos_facets
                }

    dk = tuple([ dCut(lsetp1, dt, skeleton=True, definedonelements=ba_facets[dt], deformation=lsetadap.deform)
               for dt in [NEG,POS] ])
    dGeom = tuple([ dx(definedonelements=ci.GetElementsOfType(dt)) 
                     for dt in [HASNEG,HASPOS]])
    
    # define bilinear form
    a = BilinearForm(Vh_Gamma, symmetric=True)

    # a(vh,zh)
    a += sum(mu[i] * gradv[i] * gradz[i] * dC[i] for i in [0, 1])
    a += sum( -k[i]**2 * v[i] * z[i] * dC[i] for i in [0, 1])
    if adjoint_consistent:
        a +=  average_flux_z * (v[0] - v[1]) * ds
     
    # J-alpha(u_h,v_h)
    if geom_stab_all_el: 
        sdfm = [  1.0 , 1.0 ]
    else:
        # as mentioned in paper it suffices to restricted the geometric stabilization to a band of elements 
        # around the interface (on which the meshtrafo is different from the identity)
        Vaux = L2(mesh, order=0, dirichlet=[],dgjumps=True)
        support_dfm  = GridFunction(Vaux)
        lsetadap_aux = LevelSetMeshAdaptation(mesh, order=max(2,order_geom), levelset=problem.lset)
        support_dfm.Set(IfPos(lsetadap_aux.deform.Norm(),1.0,0.0))
        for i in range(Vaux.ndof):
            if( support_dfm.vec[i]) > 1e-14:
                support_dfm.vec[i] = 1.0
        sdfm = [  support_dfm , support_dfm ]
    a += sum( alpha_stab   * h**(2*min(order,order_geom) ) * u[i] * v[i] * dGeom[i] for i in [0, 1])
    a += sum( gamma_geom   *  sdfm[i]  * h**(2*min(order,order_geom) ) * grad( u[i]) * grad(v[i]) * dGeom[i] for i in [0, 1])

    # J-Gamma(u_h,v_h)
    a += gamma_IF * (mubar/h) * (u[0] - u[1]) * (v[0] - v[1])  * ds 
    a += gamma_IF * h * jump_flux_u * jump_flux_v * ds 
    a += gamma_IF * mubar * h * jump_tangential_u * jump_tangential_v * ds  
    # J-CIP(u_h,v_h)
    a += sum( [ gamma_CIP_order * h * mu[i] * InnerProduct( (gradu[i] - gradu[i].Other()) * nF , (gradv[i] - gradv[i].Other()) * nF ) * dk[i] for i in [0,1] ]  )
    a += sum( [ gamma_GLS_order * h**2 * calL(u)[i] * calL(v)[i] * dC[i] for i in [0, 1] ] )

    # (uh,vh)_{omega} 
    a += sum( gamma_data * u[i] * v[i] * dCut(lsetp1, dt,definedon=mesh.Materials("omega"), deformation=lsetadap.deform)
            for i,dt in zip([0, 1],[NEG,POS]) )

    # a(uh,wh)
    a += sum(mu[i] * gradu[i] * gradw[i] * dC[i] for i in [0, 1])
    a += sum( -k[i]**2 * u[i] * w[i] * dC[i] for i in [0, 1])
    if adjoint_consistent:
        a +=  average_flux_w * (u[0] - u[1]) * ds

    # -s*(zh,wh) 
    a += sum(-mu[i] * gradz[i] * gradw[i] * dC[i] for i in [0, 1])

    # R.h.s.:
    f = LinearForm(Vh_Gamma)
    f += sum(coef_f[i] * w[i] * dC[i] for i in [0, 1])
    f += sum(gamma_data * solution[i] * v[i] * dCut(lsetp1, dt,definedon=mesh.Materials("omega"), deformation=lsetadap.deform)
            for i,dt in zip([0, 1],[NEG,POS]) ) 
    f += sum( [ gamma_GLS_order * h**2 * coef_f[i] * calL(v)[i] * dC[i] for i in [0, 1] ] )

    # adding noise 
    if theta_perturb != None: 
        Vh_noise_neg = Compress(Vh, GetDofsOfElements(Vh, hasneg)) 
        Vh_noise_pos = Compress(Vh, GetDofsOfElements(Vh, haspos)) 
        q_delta_neg = GridFunction(Vh_noise_neg) 
        q_delta_pos = GridFunction(Vh_noise_pos) 
        f_delta_neg = GridFunction(Vh_noise_neg)
        f_delta_pos = GridFunction(Vh_noise_pos)
        f_delta = [f_delta_neg, f_delta_pos]
        q_delta = [ q_delta_neg, q_delta_pos]
        # fill with random numbers and normalize 
        
        for i,dt in zip([0,1],[NEG,POS] ) : 
            q_delta[i].vec.FV().NumPy()[:] = np.random.rand(len( q_delta[i].vec.FV().NumPy()) )[:] 
            q_delta_norm = sqrt(Integrate(q_delta[i]**2 * dCut(lsetp1, dt ,definedon=mesh.Materials("omega")),  mesh))
            if abs(q_delta_norm) > 1e-12:
                q_delta[i].vec.FV().NumPy()[:] *= 1.0/ q_delta_norm 
        
        for i,dt in zip([0,1],[NEG,POS] ) : 
            f_delta[i].vec.FV().NumPy()[:] = np.random.rand(len( f_delta[i].vec.FV().NumPy()) )[:]
            f_norm_form =  f_delta[i]**2  * dCut(lsetp1, dt, deformation=lsetadap.deform, order=2*order)
            f_delta_norm = sqrt(Integrate(f_norm_form, mesh))
            f_delta[i].vec.FV().NumPy()[:] *= 1.0/ f_delta_norm

        if not delta_p:
            delta_p = [1,1,1]

        f += sum( delta_p[order-1] * h**(order-theta_perturb) * f_delta[i] * w[i] * dC[i] for i in [0, 1])
        f += sum( delta_p[order-1] * h**(order-theta_perturb) * q_delta  * v[i] * dCut(lsetp1, dt,definedon=mesh.Materials("omega"), deformation=lsetadap.deform)
                for i,dt in zip([0,1],[NEG,POS]) )
        f += sum( [  h**2 * h**(order-theta_perturb) *  f_delta[i] * calL(v)[i] * dC[i] for i in [0, 1] ] )

    # sets boundary data if problem well-posed (has no effect when dirichlet=[])
    with lsetadap:
        gfu.components[1].Set(solution[1], BND)

    # Setting up matrix and vector
    print("assembling linear system")
    with TaskManager():
        a.Assemble()
        f.Assemble()

    # Homogenization of boundary data and solution of linear system
    f.vec.data -= a.mat * gfu.vec
    print("Solving linera system")
    gfu.vec.data += a.mat.Inverse(Vh_Gamma.FreeDofs(),inverse=solver  )* f.vec
     
    # Export to vtk if required
    if vtk_output: 
        vtk_str = problem.lset_name + "-" + problem.problem_type + "-p{0}".format(order)+"-q{0}".format(order_geom)+"-mus({0},{1})".format(int(mu[0]),int(mu[1]))+"-ks({0},{1})".format(int(k[0]),int(k[1]))+"-lvl{0}".format(n_refs)
        #VTKOutput(ma=mesh, coefs=[lsetp1, IfPos(-lsetp1, solution[0], solution[1]) ,  IfPos(-lsetp1, sqrt((gfu.components[0]-solution[0])**2) , sqrt((gfu.components[1]-solution[1])**2)) ],
        #              names=["P1-levelset", "u", "error"],
        #              filename=vtk_str, subdivision=2).Do()

        VTKOutput(ma=mesh, coefs=[lsetp1, IfPos(-lsetp1, solution[0], solution[1]) ,  IfPos(-lsetp1, sqrt((gfu.components[0]-solution[0])**2) , sqrt((gfu.components[1]-solution[1])**2)), IfPos(-lsetp1,gfu.components[0] , gfu.components[1]) ],
                      names=["P1-levelset", "u", "error","uh"],
                      filename=vtk_str, subdivision=0).Do()
        if problem.lset_name == "ball-4-norm-convex-3D":
            domain_values = {'only_B': 0.0 , 'omega':1.0,'rest':0.0}
            values_list = [domain_values[mat]
                       for mat in mesh.GetMaterials()]
            cf_subdom = CoefficientFunction(values_list)
            vtk_str_mesh = problem.lset_name + "-" "mesh"
            VTKOutput(ma=mesh, coefs=[cf_subdom ],
                          names=["indc"],
                          filename=vtk_str_mesh, subdivision=0).Do()
            vtk_str_mesh = problem.lset_name + "-" "subdom-subdiv"
            domain_values = {'only_B': IfPos(-problem.lset,1.0, 0.0) , 'omega':1.0,'rest':0.0}
            values_list = [domain_values[mat]
                       for mat in mesh.GetMaterials()]
            cf_subdom = CoefficientFunction(values_list)
            VTKOutput(ma=mesh, coefs=[IfPos(-problem.lset,1.0, 0.0) ],
                          names=["OmegaOne"],
                          filename=vtk_str_mesh, subdivision=4).Do()


        
        if problem.lset_name == "ball-2-norm-non-convex-3D":
            print("mesh.GetMaterials() = ", mesh.GetMaterials())
            #domain_values = {'only_B': 0.0 , 'omega':1.0,'rest':0.0}
            
            domain_values = {'top': 0.0 , 'right_lid': 0.0, 'omega': 1.0, 'only_B': 0.0 }
            values_list = [domain_values[mat]
                       for mat in mesh.GetMaterials()]
            cf_subdom = CoefficientFunction(values_list)
            
            domain_values_B = {'top': 0.0 , 'right_lid': 0.0, 'omega': 0.0, 'only_B': 1.0 }
            values_list_B = [domain_values_B[mat]
                       for mat in mesh.GetMaterials()]
            cf_subdom_B = CoefficientFunction(values_list_B)


            vtk_str_mesh = problem.lset_name + "-" "mesh"
            VTKOutput(ma=mesh, coefs=[cf_subdom,cf_subdom_B  ],
                          names=["omega","onlyB"],
                          filename=vtk_str_mesh, subdivision=0).Do()
            #Draw(mesh)

        if problem.lset_name == "ball-2-norm-concentric":
            print("mesh.GetMaterials() = ", mesh.GetMaterials())
            #domain_values = {'only_B': 0.0 , 'omega':1.0,'rest':0.0}
            
            domain_values = {'rest': 0.0 , 'omega': 1.0, 'only_B': 0.0 }
            values_list = [domain_values[mat]
                       for mat in mesh.GetMaterials()]
            cf_subdom = CoefficientFunction(values_list)
            
            domain_values_B = {'rest': 0.0 , 'omega': 0.0, 'only_B': 1.0 }
            values_list_B = [domain_values_B[mat]
                       for mat in mesh.GetMaterials()]
            cf_subdom_B = CoefficientFunction(values_list_B)

            vtk_str_mesh = problem.lset_name + "-" "mesh"
            VTKOutput(ma=mesh, coefs=[cf_subdom,cf_subdom_B  ],
                          names=["omega","onlyB"],
                          filename=vtk_str_mesh, subdivision=0).Do()
            #Draw(mesh)
            #input("")

    # measure error and return
    #mat_dom = "only_B|omega"
    mat_dom = "omega|only_B"
    sol_sqr = sum([(solution[i])**2 * dCut(lsetp1, dt,definedon=mesh.Materials(mat_dom), deformation=lsetadap.deform, order=2*order) for i,dt in zip([0, 1],[NEG,POS] ) ])
    sol_norm = sqrt(Integrate(sol_sqr, mesh))
    grad_sqr = sum([ InnerProduct(sol_gradient[i], sol_gradient[i]) * dCut(lsetp1, dt,definedon=mesh.Materials(mat_dom), deformation=lsetadap.deform, order=2*order) for i,dt in zip([0, 1],[NEG,POS] ) ])
    grad_norm = sqrt(Integrate(grad_sqr, mesh))

    err_sqr = (gfuh[0]   - solution[0])**2 * dCut(lsetp1, NEG,definedon=mesh.Materials(mat_dom), deformation=lsetadap.deform, order=2*order)
    err_sqr += (gfuh[1]  - solution[1])**2 * dCut(lsetp1, POS ,definedon=mesh.Materials(mat_dom), deformation=lsetadap.deform, order=2*order)
    l2_err = sqrt(Integrate(err_sqr, mesh))
    rel_err = l2_err / sol_norm

    grad_err_sqr = InnerProduct(grad(gfuh[0]) - sol_gradient[0], grad(gfuh[0]) - sol_gradient[0] ) * dCut(lsetp1, NEG,definedon=mesh.Materials(mat_dom), deformation=lsetadap.deform, order=2*order)
    grad_err_sqr += InnerProduct(grad(gfuh[1]) - sol_gradient[1], grad(gfuh[1]) - sol_gradient[1]) * dCut(lsetp1, POS,definedon=mesh.Materials(mat_dom), deformation=lsetadap.deform, order=2*order)
    grad_err = sqrt(Integrate(grad_err_sqr, mesh))
    rel_h1_sem_err =  grad_err / grad_norm 
  
    print("L2 error : {0}, H1-semi error = {1} ".format(rel_err, rel_h1_sem_err ))
    result = {"rel-l2-err": rel_err,
              "rel-h1sem-err": rel_h1_sem_err,  
              "ndof": Vh_Gamma.ndof, 
             }

    if problem.eval_pts != []:
        uh_vals = [ ]
        u_vals = [ ]
        
        for pm in problem.eval_pts:
            if lsetp1(mesh(*pm)) < 0: 
                sub_dom_idx = 0
            else:
                sub_dom_idx = 1
            uh_vals.append( gfuh[sub_dom_idx](mesh(*pm)))
            u_vals.append( solution[sub_dom_idx](mesh(*pm)))

        result["uh-vals"] = uh_vals
        result["u-vals"] = u_vals

    return result 


def TransformToMaterialCF(mesh,fun):

    fun_vals = {'rest': fun[1],  'B_outer': fun[1], 'B_inner': fun[0], 'omega':fun[0] }
    fun_material_CF = mesh.MaterialCF(fun_vals)
    return fun_material_CF 


def SolveFitted(problem, order = 1, n_refs = 0, order_geom = 1, theta_perturb = None, order_dual = 1, stabi_dict = None, geom_stab_all_el = True, vtk_output = False,adjoint_consistent=True, delta_p=None,solver="sparsecholesky"): 

    # retrieving stabilization parameters
    if stabi_dict != None:
        gamma_CIP_order = stabi_dict["gamma-CIP"]
        gamma_GLS_order = stabi_dict["gamma-GLS"]
        alpha_stab = stabi_dict["alpha-stab"]
        gamma_data = stabi_dict["gamma-data"]
   
    # retrieve problem parameters
    
    mesh = problem.meshes[n_refs]
    mesh.Curve(order_geom)
    
    #print("mesh.GetMaterials() = ", mesh.GetMaterials())
    #print("mesh.GetBoundaries()) = ", mesh.GetBoundaries())

    mu = TransformToMaterialCF(mesh,problem.mu)
    k = TransformToMaterialCF(mesh,problem.k)
    solution = TransformToMaterialCF(mesh,problem.solution)
    sol_gradient = TransformToMaterialCF(mesh,problem.gradient)
    coef_f = TransformToMaterialCF(mesh,problem.coef_f)

    
    # background FE space for primal variable
    if problem.well_posed:
        Vh = H1(mesh, order=order, dirichlet="bc_Omega",dgjumps=True)
    else:
        Vh = H1(mesh, order=order, dirichlet=[],dgjumps=True)
    
    # continuous space for dual variable
    Vh0 = H1(mesh, order=order_dual, dirichlet="bc_Omega",dgjumps=False)
    Vh_Gamma = Vh * Vh0

    # for storing solution
    gfu = GridFunction(Vh_Gamma)
    gfuh = gfu.components
    
    # test and trial functions
    uh, zh =  Vh_Gamma.TrialFunction()
    vh, wh =  Vh_Gamma.TestFunction()
    
    gradu, gradz, gradv, gradw = [grad(fun)  for fun in [uh, zh, vh, wh]]
    gfuh = gfu.components

    # Coefficients / parameters:
    h = specialcf.mesh_size
    nF = specialcf.normal(mesh.dim)
    dF = dx(skeleton=True) 

    # differential operator for GLS
    def calL(fun):
        hesse = fun.Operator("hesse") 
        if  mesh.dim == 2:
            return -mu*hesse[0,0]-mu*hesse[1,1] - k**2*fun 
        else:
            return -mu*hesse[0,0]-mu*hesse[1,1]-mu*hesse[2,2] - k**2*fun
 
    # define bilinear form
    a = BilinearForm(Vh_Gamma, symmetric=True)

    # a(vh,zh)
    a += mu * gradv * gradz * dx
    a +=  (-1)*k**2 * vh * zh * dx 
    
    # J-CIP(u_h,v_h)
    a +=  gamma_CIP_order * h  * InnerProduct( (gradu - gradu.Other()) * nF , (gradv - gradv.Other()) * nF ) * dF
    # J-GLS
    a +=  gamma_GLS_order * h**2 * calL(uh) * calL(vh) * dx
    # Tikhonov term
    a += alpha_stab * h**(2*order) * uh * vh * dx

    # (uh,vh)_{omega} 
    a += gamma_data * uh * vh * dx(definedon=mesh.Materials("omega"))

    # a(uh,wh)
    a += mu * gradu * gradw * dx 
    a +=  (-1)*k**2 * uh * wh * dx 

    # -s*(zh,wh) 
    a += -mu * gradz * gradw * dx

    # R.h.s.:
    f = LinearForm(Vh_Gamma)
    f += coef_f * wh * dx 
    f += gamma_data * solution * vh * dx(definedon=mesh.Materials("omega"))
    f += gamma_GLS_order * h**2 * coef_f * calL(vh) * dx

    # sets boundary data if problem well-posed (has no effect when dirichlet=[])
    if problem.well_posed:
        gfu.components[0].Set(solution, BND)

    # Setting up matrix and vector
    print("assembling linear system")
    with TaskManager():
        a.Assemble()
        f.Assemble()

    # Homogenization of boundary data and solution of linear system
    f.vec.data -= a.mat * gfu.vec
    print("Solving linera system")
    gfu.vec.data += a.mat.Inverse(Vh_Gamma.FreeDofs(),inverse=solver  )* f.vec
     
    # Export to vtk if required
    #if vtk_output: 
    if False: 
        vtk_str = problem.lset_name + "-" + problem.problem_type + "-p{0}".format(order)+"-q{0}".format(order_geom)+"-mus({0},{1})".format(int(mu[0]),int(mu[1]))+"-ks({0},{1})".format(int(k[0]),int(k[1]))+"-lvl{0}".format(n_refs)
        #VTKOutput(ma=mesh, coefs=[lsetp1, IfPos(-lsetp1, solution[0], solution[1]) ,  IfPos(-lsetp1, sqrt((gfu.components[0]-solution[0])**2) , sqrt((gfu.components[1]-solution[1])**2)) ],
        #              names=["P1-levelset", "u", "error"],
        #              filename=vtk_str, subdivision=2).Do()

        VTKOutput(ma=mesh, coefs=[lsetp1, IfPos(-lsetp1, solution[0], solution[1]) ,  IfPos(-lsetp1, sqrt((gfu.components[0]-solution[0])**2) , sqrt((gfu.components[1]-solution[1])**2)), IfPos(-lsetp1,gfu.components[0] , gfu.components[1]) ],
                      names=["P1-levelset", "u", "error","uh"],
                      filename=vtk_str, subdivision=0).Do()

        #if problem.lset_name == "ball-2-norm-concentric":
        if False:
            print("mesh.GetMaterials() = ", mesh.GetMaterials())
            #domain_values = {'only_B': 0.0 , 'omega':1.0,'rest':0.0}
            
            domain_values = {'rest': 0.0 , 'omega': 1.0, 'only_B': 0.0 }
            values_list = [domain_values[mat]
                       for mat in mesh.GetMaterials()]
            cf_subdom = CoefficientFunction(values_list)
            
            domain_values_B = {'rest': 0.0 , 'omega': 0.0, 'only_B': 1.0 }
            values_list_B = [domain_values_B[mat]
                       for mat in mesh.GetMaterials()]
            cf_subdom_B = CoefficientFunction(values_list_B)

            vtk_str_mesh = problem.lset_name + "-" "mesh"
            VTKOutput(ma=mesh, coefs=[cf_subdom,cf_subdom_B  ],
                          names=["omega","onlyB"],
                          filename=vtk_str_mesh, subdivision=0).Do()
    #Draw(mesh)
    #Draw(solution,mesh,"sol")
    #Draw(gfuh[0],mesh,"uh")
    #input("")

    # measure error and return
    mat_dom = "omega|B_inner|B_outer"
    sol_sqr = (solution)**2 * dx(definedon=mesh.Materials(mat_dom))
    sol_norm = sqrt(Integrate(sol_sqr, mesh))
    grad_sqr = InnerProduct(sol_gradient, sol_gradient) * dx(definedon=mesh.Materials(mat_dom)) 
    grad_norm = sqrt(Integrate(grad_sqr, mesh))

    err_sqr = (gfuh[0]- solution)**2 * dx(definedon=mesh.Materials(mat_dom))
    l2_err = sqrt(Integrate(err_sqr, mesh))
    rel_err = l2_err / sol_norm

    grad_err_sqr = InnerProduct(grad(gfuh[0]) - sol_gradient, grad(gfuh[0]) - sol_gradient ) * dx(definedon=mesh.Materials(mat_dom))
    grad_err = sqrt(Integrate(grad_err_sqr, mesh))
    rel_h1_sem_err =  grad_err / grad_norm 
  
    print("L2 error : {0}, H1-semi error = {1} ".format(rel_err, rel_h1_sem_err ))

    result = {"rel-l2-err": rel_err,
              "rel-h1sem-err": rel_h1_sem_err,  
              "ndof": Vh_Gamma.ndof, 
             }
    
    if problem.eval_pts != []:
        uh_vals = [gfuh[0](mesh(*pm)) for pm in problem.eval_pts ]
        u_vals = [solution(mesh(*pm)) for pm in problem.eval_pts ]
        result["uh-vals"] = uh_vals
        result["u-vals"] = u_vals

    return result




def TransformUsingIndicator(mesh,fun,lset):    
    return IfPos(lset, fun[1], fun[0])


def SolveEngineering(problem, order = 1, n_refs = 0, order_geom = 1, theta_perturb = None, order_dual = 1, stabi_dict = None, geom_stab_all_el = True, vtk_output = False,adjoint_consistent=True,
                delta_p=None,solver="sparsecholesky" ): 

    # retrieving stabilization parameters
    if stabi_dict != None:
        gamma_CIP_order = stabi_dict["gamma-CIP"]
        gamma_GLS_order = stabi_dict["gamma-GLS"]
        alpha_stab = stabi_dict["alpha-stab"]
        gamma_data = stabi_dict["gamma-data"]
   
    # retrieve problem parameters
    
    mesh = problem.meshes[n_refs]
    mesh.Curve(order_geom)
    
    #print("mesh.GetMaterials() = ", mesh.GetMaterials())
    #print("mesh.GetBoundaries()) = ", mesh.GetBoundaries())

    mu = TransformUsingIndicator(mesh,problem.mu,problem.lset)
    k =  TransformUsingIndicator(mesh,problem.k,problem.lset)
    solution =  TransformUsingIndicator(mesh,problem.solution,problem.lset)
    sol_gradient = TransformUsingIndicator(mesh,problem.gradient,problem.lset)
    coef_f = TransformUsingIndicator(mesh,problem.coef_f,problem.lset)
    
    # background FE space for primal variable
    if problem.well_posed:
        Vh = H1(mesh, order=order, dirichlet="bc_Omega",dgjumps=True)
    else:
        Vh = H1(mesh, order=order, dirichlet=[],dgjumps=True)
    
    # continuous space for dual variable
    Vh0 = H1(mesh, order=order_dual, dirichlet="bc_Omega",dgjumps=False)
    Vh_Gamma = Vh * Vh0

    # for storing solution
    gfu = GridFunction(Vh_Gamma)
    gfuh = gfu.components
    
    # test and trial functions
    uh, zh =  Vh_Gamma.TrialFunction()
    vh, wh =  Vh_Gamma.TestFunction()
    
    gradu, gradz, gradv, gradw = [grad(fun)  for fun in [uh, zh, vh, wh]]
    gfuh = gfu.components

    # Coefficients / parameters:
    h = specialcf.mesh_size
    nF = specialcf.normal(mesh.dim)
    dF = dx(skeleton=True) 

    # differential operator for GLS
    def calL(fun):
        hesse = fun.Operator("hesse") 
        if  mesh.dim == 2:
            return -mu*hesse[0,0]-mu*hesse[1,1] - k**2*fun 
        else:
            return -mu*hesse[0,0]-mu*hesse[1,1]-mu*hesse[2,2] - k**2*fun
 
    # define bilinear form
    a = BilinearForm(Vh_Gamma, symmetric=True)

    # a(vh,zh)
    a += mu * gradv * gradz * dx
    a +=  (-1)*k**2 * vh * zh * dx 
    
    # J-CIP(u_h,v_h)
    a +=  gamma_CIP_order * h  * InnerProduct( (gradu - gradu.Other()) * nF , (gradv - gradv.Other()) * nF ) * dF
    # J-GLS
    a +=  gamma_GLS_order * h**2 * calL(uh) * calL(vh) * dx
    # Tikhonov term
    a += alpha_stab * h**(2*order) * uh * vh * dx

    # (uh,vh)_{omega} 
    a += gamma_data * uh * vh * dx(definedon=mesh.Materials("omega"))

    # a(uh,wh)
    a += mu * gradu * gradw * dx 
    a +=  (-1)*k**2 * uh * wh * dx 

    # -s*(zh,wh) 
    a += -mu * gradz * gradw * dx

    # R.h.s.:
    f = LinearForm(Vh_Gamma)
    f += coef_f * wh * dx 
    f += gamma_data * solution * vh * dx(definedon=mesh.Materials("omega"))
    f += gamma_GLS_order * h**2 * coef_f * calL(vh) * dx

    # sets boundary data if problem well-posed (has no effect when dirichlet=[])
    if problem.well_posed:
        gfu.components[0].Set(solution, BND)

    # Setting up matrix and vector
    print("assembling linear system")
    with TaskManager():
        a.Assemble()
        f.Assemble()

    # Homogenization of boundary data and solution of linear system
    f.vec.data -= a.mat * gfu.vec
    print("Solving linera system")
    gfu.vec.data += a.mat.Inverse(Vh_Gamma.FreeDofs(),inverse=solver  )* f.vec
     
    # Export to vtk if required
    #if vtk_output: 
    if False: 
        vtk_str = problem.lset_name + "-" + problem.problem_type + "-p{0}".format(order)+"-q{0}".format(order_geom)+"-mus({0},{1})".format(int(mu[0]),int(mu[1]))+"-ks({0},{1})".format(int(k[0]),int(k[1]))+"-lvl{0}".format(n_refs)
        #VTKOutput(ma=mesh, coefs=[lsetp1, IfPos(-lsetp1, solution[0], solution[1]) ,  IfPos(-lsetp1, sqrt((gfu.components[0]-solution[0])**2) , sqrt((gfu.components[1]-solution[1])**2)) ],
        #              names=["P1-levelset", "u", "error"],
        #              filename=vtk_str, subdivision=2).Do()

        VTKOutput(ma=mesh, coefs=[lsetp1, IfPos(-lsetp1, solution[0], solution[1]) ,  IfPos(-lsetp1, sqrt((gfu.components[0]-solution[0])**2) , sqrt((gfu.components[1]-solution[1])**2)), IfPos(-lsetp1,gfu.components[0] , gfu.components[1]) ],
                      names=["P1-levelset", "u", "error","uh"],
                      filename=vtk_str, subdivision=0).Do()

        #if problem.lset_name == "ball-2-norm-concentric":
        if False:
            print("mesh.GetMaterials() = ", mesh.GetMaterials())
            #domain_values = {'only_B': 0.0 , 'omega':1.0,'rest':0.0}
            
            domain_values = {'rest': 0.0 , 'omega': 1.0, 'only_B': 0.0 }
            values_list = [domain_values[mat]
                       for mat in mesh.GetMaterials()]
            cf_subdom = CoefficientFunction(values_list)
            
            domain_values_B = {'rest': 0.0 , 'omega': 0.0, 'only_B': 1.0 }
            values_list_B = [domain_values_B[mat]
                       for mat in mesh.GetMaterials()]
            cf_subdom_B = CoefficientFunction(values_list_B)

            vtk_str_mesh = problem.lset_name + "-" "mesh"
            VTKOutput(ma=mesh, coefs=[cf_subdom,cf_subdom_B  ],
                          names=["omega","onlyB"],
                          filename=vtk_str_mesh, subdivision=0).Do()
    #Draw(mesh)
    #Draw(solution,mesh,"sol")
    #Draw(gfuh[0],mesh,"uh")
    #input("")

    # measure error and return
    mat_dom = "omega|only_B"
    sol_sqr = (solution)**2 * dx(definedon=mesh.Materials(mat_dom))
    sol_norm = sqrt(Integrate(sol_sqr, mesh))
    grad_sqr = InnerProduct(sol_gradient, sol_gradient) * dx(definedon=mesh.Materials(mat_dom)) 
    grad_norm = sqrt(Integrate(grad_sqr, mesh))

    err_sqr = (gfuh[0]- solution)**2 * dx(definedon=mesh.Materials(mat_dom))
    l2_err = sqrt(Integrate(err_sqr, mesh))
    rel_err = l2_err / sol_norm

    grad_err_sqr = InnerProduct(grad(gfuh[0]) - sol_gradient, grad(gfuh[0]) - sol_gradient ) * dx(definedon=mesh.Materials(mat_dom))
    grad_err = sqrt(Integrate(grad_err_sqr, mesh))
    rel_h1_sem_err =  grad_err / grad_norm 
  
    print("L2 error : {0}, H1-semi error = {1} ".format(rel_err, rel_h1_sem_err ))
    result = {"rel-l2-err": rel_err,
              "rel-h1sem-err": rel_h1_sem_err,  
              "ndof": Vh_Gamma.ndof, 
             }
    
    if problem.eval_pts != []:
        uh_vals = [gfuh[0](mesh(*pm)) for pm in problem.eval_pts ]
        u_vals = [solution(mesh(*pm)) for pm in problem.eval_pts ]
        result["uh-vals"] = uh_vals
        result["u-vals"] = u_vals

    return result 
