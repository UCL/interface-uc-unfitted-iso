from netgen.geom2d import SplineGeometry, CSG2d, Circle, Rectangle
from ngsolve import *
from xfem import *
from xfem.lsetcurv import LevelSetMeshAdaptation
from math import pi
import numpy as np
from meshes import get_geometry

class interface_problem():
    def __init__(self,lset,solution,mu,k,lset_name="ball-4-norm",domain_type=None,dim=2):
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
        #self.coef_f = [-mu[i] * (solution[i].Diff(x).Diff(x)
        #               + solution[i].Diff(y).Diff(y))  - k[i]**2*solution[i]  for i in range(2)]
        self.coef_f = [-mu[i] * (solution[i].Diff(x).Diff(x)
                       + solution[i].Diff(y).Diff(y) + solution[i].Diff(z).Diff(z) )  - k[i]**2*solution[i]  for i in range(2)]
        self.well_posed = False
        self.problem_type = "ip" 
    def SetProblemType(self,well_posed):
        self.well_posed = well_posed 
        if self.well_posed:
            self.problem_type = "wp"
        else:
            self.problem_type = "ip"
    def SetDomainType(self,domain_type,ref_lvl=6):
        self.domain_type = domain_type 
        type_to_name = {"squares": "ball-4-norm-squares", 
                        "convex":"ball-4-norm-convex", 
                        "convex-3D":"ball-4-norm-convex-3D", 
                        "data-all-around":"ball-4-norm-data-all-around", 
                        "data-half": "ball-4-norm-data-half" 
                        }
        #type_to_name = {"squares": lset_name + "-squares", 
        #                "convex": lset_name + "-convex", 
        #                "data-all-around": lset_name + "-data-all-around", 
        #                "data-half": lset_name + "-data-half" 
        #                }
        if domain_type in type_to_name: 
            self.lset_name = type_to_name[domain_type] 
        for i in range(ref_lvl):
            mesh = Mesh(get_geometry(case_str= self.domain_type))
            for n in range(i):
                mesh.Refine() 
            self.meshes.append(mesh) 
    def Update(self,mu,k,solution):
        self.k = k
        self.mu = mu 
        self.solution = solution
        #self.gradient = [ CoefficientFunction((self.solution[i].Diff(x),self.solution[i].Diff(y))) for i in range(2) ] 
        if self.dim == 2:
            self.gradient = [ CoefficientFunction((self.solution[i].Diff(x),self.solution[i].Diff(y))) for i in range(2) ] 
        else:
            self.gradient = [ CoefficientFunction((self.solution[i].Diff(x),self.solution[i].Diff(y),self.solution[i].Diff(z))) for i in range(2) ] 
        #self.coef_f = [-self.mu[i] * (self.solution[i].Diff(x).Diff(x)
        #               + self.solution[i].Diff(y).Diff(y))  - self.k[i]**2*self.solution[i]  for i in range(2)]
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
    #solution = [1 + pi / 2 - sqrt(2.0) * cos(pi / 4 * r44), pi / 2 * r41]
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

# Helmholtz problem 

rho = x**2 + y**2 + z**2
#rho = x**2 + y**2
#def refsol_Helmholtz_2ball(mu,k):
#    c2 = -(k[1]/k[0])*(mu[1]/mu[0])*cos(k[1])/sin(k[0])
#    c1 = sin(k[1]) - c2*cos(k[0]) 
#    solution = [ c1 + c2 * cos( k[0] * rho),  sin( k[1] * rho) ]
#    return solution

def refsol_Helmholtz_2ball(mu,k):
    #k = [k[0]/(2*pi),k[1]/(2*pi)]
    k = [k[0],k[1]]
    c2 = -(k[1]/k[0])*(mu[1]/mu[0])*cos(k[1])/sin(k[0])
    c1 = sin(k[1]) - c2*cos(k[0]) 
    solution = [ c1 + c2 * cos( k[0] * rho),  sin( k[1] * rho) ]
    #solution = [ c1 + c2 * cos( k[0] * rho) , c1 + c2 * cos( k[0] * rho)    ]
    #solution = [ c1 + c2 * cos(  rho) , c1 + c2 * cos(  rho)    ]
    #solution = [  cos( x+y+z ) ,  cos( x+y+z)    ]
    return solution


#k = [0.05, 0.05]
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


# Pure diffusion problem in 3D 



# pure diffusion problem on flower
#r0 = 0.9 
#beta = 1/10 
#omega = 8 
#theta = atan2(x,y)
#levelset_flower = sqrt(r22) - (r0+beta*sin(omega*theta))
#eta = r22 - (r0+beta*sin(omega*theta))**2
#refsol_diffusion_flower = [ 1 + 1*cos(eta), eta**2 + 2 ]
#mu_flower = [1.0, 2.0]
#k_flower = [0.0, 0.0]
#diffusion_flower = interface_problem(lset=levelset_flower,
#                                    solution= refsol_diffusion_flower,
#                                    mu = mu_flower, 
#                                    k = k_flower)





