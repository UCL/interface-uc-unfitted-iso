from netgen.geom2d import SplineGeometry, CSG2d, Circle, Rectangle
from ngsolve import *
from xfem import *
from xfem.lsetcurv import LevelSetMeshAdaptation
from math import pi
from netgen.csg import *
import numpy as np

def get_geometry(case_str="nonconvex"):
    if case_str == "squares":
        geo = CSG2d()
        omega_dom = Rectangle( pmin=(-0.5,-0.5), pmax=(0.5,0.5), mat="omega", bc="bc_omega")
        B_dom = Rectangle( pmin=(-1.25,-1.25), pmax=(1.25,1.25), mat="B", bc="bc_B")
        full_dom = Rectangle( pmin=(-1.5,-1.5), pmax=(1.5,1.5), mat="full", bc="bc_Omega")

        only_B = B_dom - omega_dom
        only_B.Mat("only_B")
        rest = full_dom - B_dom 

        # add top level objects to geometry
        geo.Add(omega_dom)
        geo.Add(only_B)
        geo.Add(rest)

        # generate mesh
        m = geo.GenerateMesh(maxh=0.25)
        return m

    if case_str == "convex":
        geo = SplineGeometry() 
        # data domain
        p1 = geo.AppendPoint (-1.5,-1.5)
        p2 = geo.AppendPoint (1.5,-1.5)
        p3 = geo.AppendPoint (1.5,1.25)
        p4 = geo.AppendPoint (1.25,1.25)
        p5 = geo.AppendPoint (1.25,-1.25)
        p6 = geo.AppendPoint (-1.25,-1.25)
        p7 = geo.AppendPoint (-1.25,1.25)
        p8 = geo.AppendPoint (-1.5,1.25)

        # omega
        geo.Append (["line", p1, p2], leftdomain=1, rightdomain=0,bc="bc_Omega")
        geo.Append (["line", p2, p3], leftdomain=1, rightdomain=0,bc="bc_Omega")
        geo.Append (["line", p3, p4], leftdomain=1, rightdomain=3)
        geo.Append (["line", p4, p5], leftdomain=1, rightdomain=2)
        geo.Append (["line", p5, p6], leftdomain=1, rightdomain=2)
        geo.Append (["line", p6, p7], leftdomain=1, rightdomain=2)
        geo.Append (["line", p7, p8], leftdomain=1, rightdomain=3)
        geo.Append (["line", p8, p1], leftdomain=1, rightdomain=0,bc="bc_Omega")

        # B 
        geo.Append (["line", p4, p7], leftdomain=2, rightdomain=3)

        # rest 
        p11 = geo.AppendPoint(1.5,1.5)
        p12 = geo.AppendPoint(-1.5,1.5)
        geo.Append (["line", p3, p11], leftdomain=3, rightdomain=0,bc="bc_Omega")
        geo.Append (["line", p11, p12], leftdomain=3, rightdomain=0,bc="bc_Omega")
        geo.Append (["line", p12, p8], leftdomain=3, rightdomain=0,bc="bc_Omega")

        geo.SetMaterial(1, "omega")
        geo.SetMaterial(2, "only_B")
        geo.SetMaterial(3, "rest")

        m = geo.GenerateMesh(maxh=0.25)
        return m

    if case_str == "data-all-around":
        geo = SplineGeometry() 
        # data domain
        p1 = geo.AppendPoint (-1.5,-1.5)
        p2 = geo.AppendPoint (1.5,-1.5)
        p3 = geo.AppendPoint (1.5,1.25)
        p4 = geo.AppendPoint (1.25,1.25)
        p5 = geo.AppendPoint (1.25,-1.25)
        p6 = geo.AppendPoint (-1.25,-1.25)
        p7 = geo.AppendPoint (-1.25,1.25)
        p8 = geo.AppendPoint (-1.5,1.25)
        p11 = geo.AppendPoint(1.5,1.5)
        p12 = geo.AppendPoint(-1.5,1.5)

        # omega
        geo.Append (["line", p1, p2], leftdomain=1, rightdomain=0,bc="bc_Omega")
        geo.Append (["line", p2, p11], leftdomain=1, rightdomain=0,bc="bc_Omega")
        geo.Append (["line", p11, p12], leftdomain=1, rightdomain=0,bc="bc_Omega")
        geo.Append (["line", p12, p1], leftdomain=1, rightdomain=0,bc="bc_Omega")

        # only B 
        geo.Append (["line", p6, p5], leftdomain=2, rightdomain=1)
        geo.Append (["line", p5, p4], leftdomain=2, rightdomain=1)
        geo.Append (["line", p4, p7], leftdomain=2, rightdomain=1)
        geo.Append (["line", p7, p6], leftdomain=2, rightdomain=1)

        geo.SetMaterial(1, "omega")
        geo.SetMaterial(2, "only_B")

        m = geo.GenerateMesh(maxh=0.25)
        return m

    if case_str == "data-half":
        geo = SplineGeometry() 
        p1 = geo.AppendPoint (-1.5,-1.5)
        p2 = geo.AppendPoint (1.5,-1.5)
        p3 = geo.AppendPoint (1.5,-1.25)
        p4 = geo.AppendPoint (1.25,-1.25)
        p5 = geo.AppendPoint (-1.25,-1.25)
        p6 = geo.AppendPoint (-1.5,-1.25)
        p7 = geo.AppendPoint (1.25,1.25)
        p8 = geo.AppendPoint (-1.25,1.25)
        p9 = geo.AppendPoint (1.5,1.5)
        p10 = geo.AppendPoint (-1.5,1.5)
        p11 = geo.AppendPoint (1.5,0.0)
        p12 = geo.AppendPoint (1.25,0.0)
        p13 = geo.AppendPoint (-1.25,0.0)
        p14 = geo.AppendPoint (-1.5,0.0)

        # omega
        geo.Append (["line", p1, p2], leftdomain=1, rightdomain=0,bc="bc_Omega")
        geo.Append (["line", p2, p11], leftdomain=1, rightdomain=0,bc="bc_Omega") 
        geo.Append (["line", p11, p12], leftdomain=1, rightdomain=3)
        geo.Append (["line", p12, p4], leftdomain=1, rightdomain=2)
        geo.Append (["line", p4, p5], leftdomain=1, rightdomain=2)
        geo.Append (["line", p5, p13], leftdomain=1, rightdomain=2)
        geo.Append (["line", p13, p14], leftdomain=1, rightdomain=3)
        geo.Append (["line", p14, p1], leftdomain=1, rightdomain=0,bc="bc_Omega") 

        # only B
        geo.Append (["line", p12, p7], leftdomain=2, rightdomain=3)
        geo.Append (["line", p7, p8], leftdomain=2, rightdomain=3)
        geo.Append (["line", p8, p13], leftdomain=2, rightdomain=3)

        # rest 
        geo.Append (["line", p11, p9], leftdomain=3, rightdomain=0,bc="bc_Omega")
        geo.Append (["line", p9, p10], leftdomain=3, rightdomain=0,bc="bc_Omega") 
        geo.Append (["line", p10, p14], leftdomain=3, rightdomain=0,bc="bc_Omega")

        geo.SetMaterial(1, "omega")
        geo.SetMaterial(2, "only_B")
        geo.SetMaterial(3, "rest")

        m = geo.GenerateMesh(maxh=0.25)
        return m

    if case_str == "convex-3D":
        
        print("Creating 3D mesh")
        geo = CSGeometry()
        cube = OrthoBrick( Pnt(-1.75,-1.75,-1.75), Pnt(1.75,1.75,1.2) )
        top = OrthoBrick( Pnt(-1.75,-1.75,1.2), Pnt(1.75,1.75,1.75) )
        top.mat("rest") 
        
        cube.bc('bc_Omega')
        inside = OrthoBrick( Pnt(-1.2,-1.2,-1.2), Pnt(1.2,1.2,1.2) )
        inside.mat("only_B")
        omega = cube - inside 
        omega.mat("omega")
        geo.Add(omega)
        geo.Add(inside)
        geo.Add(top)
        
        m = geo.GenerateMesh(maxh=1.15)
        return m
