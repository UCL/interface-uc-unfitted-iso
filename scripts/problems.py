from netgen.geom2d import SplineGeometry, CSG2d, Circle, Rectangle
from ngsolve import *
from xfem import *
from xfem.lsetcurv import LevelSetMeshAdaptation
from math import pi
import numpy as np
from meshes import get_geometry

class interface_problem():
    def __init__(self,lset,solution,mu,k,lset_name="ball-4-norm",domain_type=None,dim=2, pre_refine_lset = False, hlist = None):
        self.lset = lset 
        self.lset_name = lset_name 
        self.solution = solution 
        self.dim = dim
        if self.dim == 2:
            self.gradient = [ CoefficientFunction((self.solution[i].Diff(x),self.solution[i].Diff(y))) for i in range(2) ] 
        else:
            self.gradient = [ CoefficientFunction((self.solution[i].Diff(x),self.solution[i].Diff(y),self.solution[i].Diff(z))) for i in range(2) ] 
        self.mu = mu
        self.k = k 
        self.meshes = [ ]
        self.domain_type = domain_type 
        self.coef_f = [-mu[i] * (solution[i].Diff(x).Diff(x)
                       + solution[i].Diff(y).Diff(y) + solution[i].Diff(z).Diff(z) )  - k[i]**2*solution[i]  for i in range(2)]
        self.well_posed = False
        self.problem_type = "ip"
        self.eval_pts = []
        self.pre_refine_lset = pre_refine_lset
        self.hlist = hlist
    def SetProblemType(self,well_posed):
        self.well_posed = well_posed 
        if self.well_posed:
            self.problem_type = "wp"
        else:
            self.problem_type = "ip"
    def SetDomainType(self,domain_type,ref_lvl=6):
        self.domain_type = domain_type 
        type_to_name = {"squares": "ball-4-norm-squares", 
                        "squares-easy": "ball-4-norm-squares", 
                        "convex":"ball-4-norm-convex", 
                        "convex-3D":"ball-4-norm-convex-3D", 
                        "data-all-around":"ball-4-norm-data-all-around", 
                        "data-half": "ball-4-norm-data-half", 
                        "non-convex-3D": "ball-2-norm-non-convex-3D",
                        "concentric-3D": "ball-2-norm-concentric",
                        "concentric-3D-fitted": "ball-2-concentric-fitted",
                        "convex-3D-concentric": "ball-2-concentric-convex",
                        "convex-3D-fitted": "ball-2-convex-fitted"
                        }
        if domain_type in type_to_name: 
            self.lset_name = type_to_name[domain_type] 
            if not self.hlist:
                for i in range(ref_lvl):
                    mesh = Mesh(get_geometry(case_str= self.domain_type))
                    #print("mesh.GetMaterials() = ", mesh.GetMaterials())
                    #print("mesh.GetBoundaries()) = ", mesh.GetBoundaries())
                    #Draw(mesh)
                    #input("")
                    if self.pre_refine_lset:
                        for j in range(1):
                            marks = [ ]
                            # refinement criterion: check if centroid is close to zero level set  
                            #centroids = [] 
                            for el in mesh.Elements():
                                vs = el.vertices
                                L = len(el.vertices)
                                c_x = 0
                                c_y = 0
                                c_z = 0
                                for vert in el.vertices:
                                    vx,vy,vz = mesh.ngmesh.Points()[vert.nr+1].p
                                    c_x = c_x + vx
                                    c_y = c_y + vy
                                    c_z = c_z + vz
                                c_x = c_x / L
                                c_y = c_y / L
                                c_z = c_z / L
                                #centroids.append( (c_x,c_y,c_z) )
                                centroid = (c_x,c_y,c_z)
                                lset_val = self.lset(mesh(*centroid))
                                #print("Centroid = {0}, lset value = {1}".format( centroid,  lset_val ))
                                if  abs(lset_val) <  self.pre_refine_lset:
                                    marks.append(True)
                                else:
                                    marks.append(False)
                            #print("Number of elements to refine = ",sum(marks))
                            for el in mesh.Elements():        
                                mesh.SetRefinementFlag(el,marks[el.nr])
                            mesh.Refine()
                    for n in range(i):
                        mesh.Refine() 
                    self.meshes.append(mesh)
            else:
                for maxh in self.hlist:
                    self.meshes.append(Mesh(get_geometry(case_str= self.domain_type,maxh=maxh)))
    def Update(self,mu,k,solution):
        self.k = k
        self.mu = mu 
        self.solution = solution
        if self.dim == 2:
            self.gradient = [ CoefficientFunction((self.solution[i].Diff(x),self.solution[i].Diff(y))) for i in range(2) ] 
        else:
            self.gradient = [ CoefficientFunction((self.solution[i].Diff(x),self.solution[i].Diff(y),self.solution[i].Diff(z))) for i in range(2) ] 
        self.coef_f = [-self.mu[i] * (self.solution[i].Diff(x).Diff(x)
                       + self.solution[i].Diff(y).Diff(y) + self.solution[i].Diff(z).Diff(z) )  - self.k[i]**2*self.solution[i]  for i in range(2)]


# pure diffusion problem on 4ball 
r22 = x**2 + y**2
r44 = x**4 + y**4
r66 = x**6 + y**6
r41 = sqrt(sqrt(r44))
r4m3 = 1.0 / r41**3

def refsol_diffusion(mu):
    c1 = (1/sqrt(2.0))*(1+pi*mu[0]/mu[1]) 
    c2 = 1
    c3 = (mu[0]/mu[1])*(pi/sqrt(2.0))
    solution = [c1 - c2 * cos(pi / 4 * r44), c3 * r41]
    return solution

levelset = r41 - 1.0
mu = [1.0, 2.0]
k = [0.0, 0.0]
diffusion_4ball = interface_problem(lset=levelset,
                                    solution= refsol_diffusion(mu),
                                    mu = mu, 
                                    k = k)

#helmholtz_4ball 
def refsol_Helmholtz(mu,k):
    c3 = 1
    c2 = -c3 * (k[1]/k[0])*(mu[1]/mu[0]) * (cos(k[1])/sin(k[0]))
    c1 = c3*sin(k[1]) - c2*cos(k[0])
    solution = [ c1 + c2 * cos( k[0] * r44), c3 *  sin( k[1] * r44) ]
    return solution
k = [4.0, 4.0]
helmholtz_4ball = interface_problem(lset=levelset,
                                    solution=refsol_Helmholtz(mu,k),
                                    mu = mu, 
                                    k = k)

rho = x**2 + y**2 + z**2

def refsol_Helmholtz_2ball(mu,k):
    #k = [k[0]/(2*pi),k[1]/(2*pi)]
    k = [k[0],k[1]]
    c2 = -(k[1]/k[0])*(mu[1]/mu[0])*cos(k[1])/sin(k[0])
    c1 = sin(k[1]) - c2*cos(k[0]) 
    solution = [ c1 + c2 * cos( k[0] * rho),  sin( k[1] * rho) ]
    return solution


def refsol_Helmholtz_2ball(mu,k):
    #alpha = 1/0.05
    alpha = 1
    #k = [k[0]/(2*pi),k[1]/(2*pi)]
    k = [k[0],k[1]]
    #c2 = -(k[1]/k[0])*(mu[1]/mu[0])*cos(k[1])/sin(k[0])
    #c1 = sin(k[1]) - c2*cos(k[0]) 
    
    c1 = -k[1]*mu[1]*cos(k[1]) / ( k[0]*mu[0]*sin(k[0]) )
    c2 =  sin(k[1]) - c1 * cos(k[0])
    rr = sqrt(rho)
    phi_rho = rho-1.0
    chi = exp(-phi_rho**2*alpha) 
    solution = [ (c1*cos(k[0]*rr) + c2) * chi , sin(k[1]*rr) * chi  ]
    

    #bump = exp( -(x-0.0)**2/0.025 - (y-0.0)**2/0.025  - (z-1.0)**2/0.025 ) 
    #solution = [ bump ,  bump ]
    
    #solution = [ exp(-phi_rho**2/0.05) , exp(-phi_rho**2/0.05)   ]
    #solution = [ cos(k[0]*rho)*exp(-phi_rho**2/0.05) , sin(k[1]*rho)*exp(-phi_rho**2/0.05)   ]
    return solution


k = [4.0,4.0]
levelset_2ball = sqrt(rho) - 1.0 
helmholtz_2ball = interface_problem(lset=levelset_2ball,
                                    solution=refsol_Helmholtz_2ball(mu,k),
                                    mu = mu, 
                                    k = k,
                                    dim=3)

r44_3D = x**4 + y**4 + z**4
r41_3D = sqrt(sqrt(r44_3D))
levelset_3D = r41_3D - 1.0 
def refsol_diffusion_3D(mu):
    c1 = (1/sqrt(2.0))*(1+pi*mu[0]/mu[1]) 
    c2 = 1
    c3 = (mu[0]/mu[1])*(pi/sqrt(2.0))
    solution = [c1 - c2 * cos(pi / 4 * r44_3D), c3 * r41_3D]
    return solution

mu = [1,2]
k = [0,0]
diffusion_3D = interface_problem(lset = levelset_3D,
                                    solution = refsol_diffusion_3D(mu),
                                    mu = mu, 
                                    k = k,
                                    dim=3)


mu = [1,1]
k = [1,1]
helmholtz_3D_ball = interface_problem(lset = levelset_2ball,
                                    solution = refsol_Helmholtz_2ball(mu,k),
                                    mu = mu, 
                                    k = k,
                                    dim=3)




