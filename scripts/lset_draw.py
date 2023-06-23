# this file contains some auxiliary functions for drawing meshes
# (it is not needed to reproduce any of the numerical results) 
from netgen.geom2d import SplineGeometry, CSG2d, Circle, Rectangle
from ngsolve import *
from xfem import *
from xfem.lsetcurv import LevelSetMeshAdaptation
from math import pi
from scipy import optimize
import numpy as np
#import matplotlib
#matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from meshes import get_geometry


r22 = x**2 + y**2
r44 = x**4 + y**4
r66 = x**6 + y**6
r41 = sqrt(sqrt(r44))

levelset = r41 - 1.0


ddx = 10


def draw_mesh_tikz(fes,mesh,name,case_str):
    ddx = 10
    rainbow = ["cyan","white"]
    file = open("{0}.tex".format(name),"w+")
    file.write("\\documentclass{standalone} \n")
    file.write("\\usepackage{xr} \n")
    file.write("\\usepackage{tikz} \n")
    file.write("\\usepackage{pgfplots} \n")
    file.write("\\usepackage{xcolor} \n")
    file.write("\\usepackage{} \n")   
    file.write("\\usetikzlibrary{shapes,arrows,shadows,snakes,calendar,matrix,spy,backgrounds,folding,calc,positioning,patterns,patterns.meta} \n")
    file.write("\\selectcolormodel{cmyk}  \n")  
    file.write("\\definecolor{orange}{cmyk}{0,0.45,0.91,0} \n")  
    file.write("\\definecolor{brightblue}{cmyk}{0.92,0,0.15,0.05} \n")  
    file.write("\\definecolor{richred}{cmyk}{0,1,0.62,0} \n")  
    file.write("\\definecolor{yellow}{cmyk}{0,0.25,0.95,0} \n")  
    file.write("\\begin{document} \n")
    file.write("\\begin{tikzpicture}[scale = 1.0] \n")  
    file.write("\\pgfresetboundingbox \n")
    file.write("\\path[use as bounding box,draw,black,ultra thick] (-{0},-{0}) rectangle ({0},{0}); \n".format(ddx*1.5) )
    for el in fes.Elements():
        x_coords = [] 
        y_coords = []
        coords = []
        for vert in el.vertices:
            vx,vy,vz = mesh.ngmesh.Points()[vert.nr+1].p
            coords.append((ddx*vx,ddx*vy))
            x_coords.append(vx)
            y_coords.append(vy)
        x_coords = np.array(x_coords)
        y_coords = np.array(y_coords)
        if len(coords) == 3:
            if case_str == "data-all-around-omega":
                if np.all( x_coords <= -1.25) or np.all( x_coords >= 1.25) or  np.all( y_coords <= -1.25) or np.all( y_coords >= 1.25):
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("yellow",coords[0],coords[1],coords[2] ))
                else:
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("white",coords[0],coords[1],coords[2] ))
            if case_str == "data-all-around-B":
                if True:
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("brightblue",coords[0],coords[1],coords[2] ))
                else:
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("white",coords[0],coords[1],coords[2] ))
            if case_str == "convex-omega":
                if  np.all( x_coords <= -1.25) and (  np.all( y_coords <= 1.25)): 
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("yellow",coords[0],coords[1],coords[2] ))
                elif np.all( x_coords >= 1.25) and (  np.all( y_coords <= 1.25)): 
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("yellow",coords[0],coords[1],coords[2] ))
                elif  (  np.all( y_coords <= -1.25)): 
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("yellow",coords[0],coords[1],coords[2] ))
                else:
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("white",coords[0],coords[1],coords[2] ))
            if case_str == "convex-B":
                if (  np.all( y_coords <= 1.25)):
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("brightblue",coords[0],coords[1],coords[2] ))
                else:
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("white",coords[0],coords[1],coords[2] ))
            if case_str == "data-half-omega":
                if  np.all( x_coords <= -1.25) and (  np.all( y_coords <= 0.0)): 
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("yellow",coords[0],coords[1],coords[2] ))
                elif np.all( x_coords >= 1.25) and (  np.all( y_coords <= 0.0)): 
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("yellow",coords[0],coords[1],coords[2] ))
                elif  (  np.all( y_coords <= -1.25)): 
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("yellow",coords[0],coords[1],coords[2] ))
                else:
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("white",coords[0],coords[1],coords[2] ))
            if case_str == "data-half-B":
                if  np.all( x_coords <= -1.25) and (  np.all( y_coords <= 0.0)): 
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("brightblue",coords[0],coords[1],coords[2] ))
                elif np.all( x_coords >= 1.25) and (  np.all( y_coords <= 0.0)): 
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("brightblue",coords[0],coords[1],coords[2] ))
                elif  (  np.all( y_coords <= -1.25)): 
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("brightblue",coords[0],coords[1],coords[2] ))
                elif np.all( x_coords >= -1.25) and np.all( x_coords <= 1.25) and  np.all( y_coords >= -1.25) and np.all( y_coords <= 1.25): 
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("brightblue",coords[0],coords[1],coords[2] ))
                else:
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("white",coords[0],coords[1],coords[2] ))
            if case_str == "squares-omega" or "squares-omega-B":
                if  np.all( x_coords >= -0.5) and np.all( x_coords <= 0.5) and  np.all( y_coords >= -0.5) and np.all( y_coords <= 0.5):  
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("yellow",coords[0],coords[1],coords[2] ))
                else:
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("white",coords[0],coords[1],coords[2] ))
            if case_str == "squares-B":
                if  np.all( x_coords >= -1.25) and np.all( x_coords <= 1.25) and  np.all( y_coords >= -1.25) and np.all( y_coords <= 1.25):  
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("brightblue",coords[0],coords[1],coords[2] ))
                else:
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("white",coords[0],coords[1],coords[2] ))

            if case_str == "convex-omega-B":
                if  np.all( x_coords <= -1.25) and (  np.all( y_coords <= 1.25)): 
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("yellow",coords[0],coords[1],coords[2] ))
                elif np.all( x_coords >= 1.25) and (  np.all( y_coords <= 1.25)): 
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("yellow",coords[0],coords[1],coords[2] ))
                elif  (  np.all( y_coords <= -1.25)): 
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=0.8] {1} -- {2} -- {3} -- cycle; \n".format("yellow",coords[0],coords[1],coords[2] ))
                else:
                    file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=1.0] {1} -- {2} -- {3} -- cycle; \n".format("white",coords[0],coords[1],coords[2] ))


        if len(coords) == 4:
            file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=1] {1} -- {2} -- {3} -- {4} -- cycle; \n".format("white",coords[0],coords[1],coords[2],coords[3]))
            

    #s1 = "\\draw[line width=2mm,draw =orange,->] (0.0, 0.0) -- node[at end,above]{ \\resizebox{ .125\\linewidth}{!}{\\itshape \\textcolor{orange}{$r$} } }"
    #s1 += "( 0.0, {0}); \n".format( 1.2*ddx)
    #file.write(s1)

    #file.write("\\draw[line width=2mm,draw =mLightBrown] (0,0) circle ({0}); \n".format( 0.326*dx))
    #s1 = "\\draw[line width=1mm,draw =black,dashed] (0.0, 0.0) -- node[midway,above]{ \\resizebox{ .125\\linewidth}{!}{\\itshape $r_{*}$} }"
    #s1 += "( {0}, 0.0); \n".format( 0.326*dx)
    #file.write(s1 )
    #\draw[pattern={dots},pattern color=orange]
    #  (0,0) rectangle +(1,1);
    if case_str == "data-all-around-omega":
        file.write("\\draw[white ] ({0},{1})  node[fill=white,above]{{  \\resizebox{{ .25\\linewidth}}{{!}}{{ \\textcolor{{black}}{{$\omega$}} }} }};   \n".format( 0.0*ddx ,-1.45*ddx  ))
    if case_str == "data-all-around-B":
        file.write("\\draw[white ] ({0},{1})  node[fill=white,above]{{  \\resizebox{{ .25\\linewidth}}{{!}}{{ \\textcolor{{black}}{{$B$}} }} }};   \n".format( 0.0*ddx ,0.0*ddx  ))
    if case_str == "convex-omega" or case_str == "data-half-omega":
        file.write("\\draw[white ] ({0},{1})  node[fill=white,above]{{  \\resizebox{{ .25\\linewidth}}{{!}}{{ \\textcolor{{black}}{{$\omega$}} }} }};   \n".format( 0.0*ddx ,-1.45*ddx  ))
    if case_str == "convex-B" or case_str == "data-half-B" or case_str == "squares-B": 
        file.write("\\draw[white ] ({0},{1})  node[fill=white,above]{{  \\resizebox{{ .25\\linewidth}}{{!}}{{ \\textcolor{{black}}{{$B$}} }} }};   \n".format( 0.0*ddx ,0.0*ddx  ))
    if case_str == "squares-omega": 
        file.write("\\draw[white ] ({0},{1})  node[fill=white,above]{{  \\resizebox{{ .25\\linewidth}}{{!}}{{ \\textcolor{{black}}{{$\omega$}} }} }};   \n".format( 0.0*ddx ,0.0*ddx  ))
    if case_str == "convex-omega-B":
        file.write("\\draw[pattern={{Dots[ distance={{30pt}},radius={{4pt}}]  }}, pattern color = magenta ] ({0},{1}) rectangle ({2},{3}); \n".format(-1.5*ddx,-1.5*ddx,1.5*ddx,1.25*ddx)) 
        file.write("\\draw[white ] ({0},{1})  node[fill=white,above]{{  \\resizebox{{ .25\\linewidth}}{{!}}{{ \\textcolor{{magenta}}{{$B$}} }} }};   \n".format( 0.0*ddx ,0.0*ddx  ))
        file.write("\\draw[white ] ({0},{1})  node[fill=white,above]{{  \\resizebox{{ .25\\linewidth}}{{!}}{{ \\textcolor{{yellow}}{{$\omega$}} }} }};   \n".format( 0.0*ddx ,-1.45*ddx  ))
    if case_str == "convex-omega-B":
        file.write("\\draw[pattern={{Dots[ distance={{30pt}},radius={{4pt}}]  }}, pattern color = magenta ] ({0},{1}) rectangle ({2},{3}); \n".format(-1.5*ddx,-1.5*ddx,1.5*ddx,1.25*ddx)) 
        file.write("\\draw[white ] ({0},{1})  node[fill=white,above]{{  \\resizebox{{ .25\\linewidth}}{{!}}{{ \\textcolor{{magenta}}{{$B$}} }} }};   \n".format( 0.0*ddx ,0.0*ddx  ))
        file.write("\\draw[white ] ({0},{1})  node[fill=white,above]{{  \\resizebox{{ .25\\linewidth}}{{!}}{{ \\textcolor{{yellow}}{{$\omega$}} }} }};   \n".format( 0.0*ddx ,-1.45*ddx  ))
    if case_str == "squares-omega-B":
        file.write("\\draw[pattern={{Dots[ distance={{30pt}},radius={{4pt}}]  }}, pattern color = magenta ] ({0},{1}) rectangle ({2},{3}); \n".format(-1.25*ddx,-1.25*ddx,1.25*ddx,1.25*ddx)) 
        file.write("\\draw[white ] ({0},{1})  node[fill=white,above]{{  \\resizebox{{ .25\\linewidth}}{{!}}{{ \\textcolor{{magenta}}{{$B$}} }} }};   \n".format( 0.0*ddx ,1.1*ddx  ))
        file.write("\\draw[white ] ({0},{1})  node[fill=white,above]{{  \\resizebox{{ .25\\linewidth}}{{!}}{{ \\textcolor{{yellow}}{{$\omega$}} }} }};   \n".format( 0.0*ddx ,0.0*ddx  ))


    if case_str == "convex-omega-B" or "squares-omega-B":
        file.write("\\draw[ ] ({0},{1})  node[fill=white,above]{{  \\resizebox{{ .25\\linewidth}}{{!}}{{ \\textcolor{{green!90!black}}{{$\Gamma$}} }} }};   \n".format( 0.0*ddx ,-0.95*ddx  ))
        file.write("\\begin{{axis}}[view={{0}}{{90}},  anchor=origin,  disabledatascaling, at={{(0pt,0pt)}}, x={0}cm,y={0}cm,z={0}cm, hide axis ] \n".format(ddx))
        file.write("\\addplot3 [line width=10pt, \n") 
        file.write(" contour lua={levels={0},labels=false, draw color=green!90!black}, samples=200 ] {x^4 + y^4 - 1}; \n")
    else:
        file.write("\\draw[ ] ({0},{1})  node[above]{{  \\resizebox{{ .25\\linewidth}}{{!}}{{ \\textcolor{{richred}}{{$\Gamma$}} }} }};   \n".format( 0.0*ddx ,-0.95*ddx  ))
        file.write("\\begin{{axis}}[view={{0}}{{90}},  anchor=origin,  disabledatascaling, at={{(0pt,0pt)}}, x={0}cm,y={0}cm,z={0}cm, hide axis ] \n".format(ddx))
        file.write("\\addplot3 [ultra thick, \n") 
        file.write(" contour lua={levels={0},labels=false, draw color=richred}, samples=200 ] {x^4 + y^4 - 1}; \n")

    
    file.write(" \end{axis} \n ")
    file.write("\\end{tikzpicture} \n") 
    file.write("\\end{document} \n")           
    file.close()

# draw convex mesh 
def draw_based_on_geom(case_str="convex"):
    mesh = Mesh(get_geometry(case_str= case_str))
    lsetadap = LevelSetMeshAdaptation(mesh, order=1, levelset=levelset)
    lsetp1 = lsetadap.lset_p1
    Vh = H1(mesh, order=1, dirichlet=[],dgjumps=True)
    Vh0 = H1(mesh, order=1, dirichlet="bc_Omega",dgjumps=False)
    ci = CutInfo(mesh, lsetp1)
    hasneg = ci.GetElementsOfType(HASNEG)
    haspos = ci.GetElementsOfType(HASPOS)
    hasif = ci.GetElementsOfType(IF)
    Vh_Gamma = Compress(Vh, GetDofsOfElements(Vh, hasneg)) \
              * Compress(Vh, GetDofsOfElements(Vh, haspos)) \
              * Vh0
    draw_mesh_tikz(Vh_Gamma,mesh, name= case_str+"-omega-B",case_str=case_str+"-omega-B")

draw_based_on_geom(case_str="convex")
draw_based_on_geom(case_str="squares")


levelset = sqrt(r22) - 1.0


def draw_mesh_simple(name="dummy",case_str="geom",draw_lset=True):
    ll, ur = (-1.5, -1.5), (1.5, 1.5)
    square = SplineGeometry()
    square.AddRectangle(ll, ur, bc=1)
    mesh = Mesh(square.GenerateMesh(maxh=0.2, quad_dominated=False))
    fes = H1(mesh, order=1, dirichlet=[])
    
    ddx = 10
    rainbow = ["cyan","white"]
    file = open("../plots/meshes/{0}.tex".format(name),"w+")
    file.write("\\documentclass{standalone} \n")
    file.write("\\usepackage{xr} \n")
    file.write("\\usepackage{tikz} \n")
    file.write("\\usepackage{pgfplots} \n")
    file.write("\\usepackage{xcolor} \n")
    file.write("\\usepackage{} \n")   
    file.write("\\usetikzlibrary{shapes,arrows,shadows,snakes,calendar,matrix,spy,backgrounds,folding,calc,positioning,patterns,hobby} \n")
    file.write("\\selectcolormodel{cmyk}  \n")  
    file.write("\\definecolor{orange}{cmyk}{0,0.45,0.91,0} \n")  
    file.write("\\definecolor{brightblue}{cmyk}{0.92,0,0.15,0.05} \n")  
    file.write("\\definecolor{richred}{cmyk}{0,1,0.62,0} \n")  
    file.write("\\definecolor{yellow}{cmyk}{0,0.25,0.95,0} \n")  
    file.write("\\begin{document} \n")
    file.write("\\begin{tikzpicture}[scale = 1.0,use Hobby shortcut] \n")  
    file.write("\\pgfresetboundingbox \n")
    file.write("\\path[use as bounding box,draw,black,ultra thick] (-{0},-{0}) rectangle ({0},{0}); \n".format(ddx*1.5) )

    for el in fes.Elements():
        x_coords = [] 
        y_coords = []
        coords = []
        for vert in el.vertices:
            vx,vy,vz = mesh.ngmesh.Points()[vert.nr+1].p
            coords.append((ddx*vx,ddx*vy))
            x_coords.append(vx)
            y_coords.append(vy)
        x_coords = np.array(x_coords)
        y_coords = np.array(y_coords)
        if len(coords) == 3:
            file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=1] {1} -- {2} -- {3} -- cycle; \n".format("white",coords[0],coords[1],coords[2] ))
        if len(coords) == 4:
            file.write("\\draw[line width=0.01mm,draw =black, fill={0},fill opacity=1] {1} -- {2} -- {3} -- {4} -- cycle; \n".format("white",coords[0],coords[1],coords[2],coords[3]))
    

    file.write("\\path  \n")
    file.write("(0.0,0.0) coordinate (z0) \n")
    file.write("(2.0,-3.0) coordinate (z1) \n")
    file.write("(6.0,-6.0) coordinate (z2)\n")
    file.write("(10,-2) coordinate (z3)\n")
    file.write("(0.0,10) coordinate (z4)\n")
    file.write("(-10,-2) coordinate (z5)\n")
    file.write("(-6,-6) coordinate (z6)\n")
    file.write("(-2,-3.0) coordinate (z7);\n")
    
    if draw_lset:
        file.write("\draw[richred, line width=2.5mm,closed] (z0) .. (z1) .. (z3) .. (z4) .. (z5) .. (z6) .. (z7); \n") 
        file.write("\\draw[] ( 0.0 , 2.0  )  node[fill=white,above]{{   \\resizebox{ .35\\linewidth}{!}{ \\textcolor{black}{$\Omega_1$}}     }}; \n")
        file.write("\\draw[] ( -13 , 4.4  )  node[fill=white,above]{{ \\resizebox{ .35\\linewidth}{!}{ \\textcolor{black}{$\Omega_2$ }}      }}; \n")
        file.write("\\draw[] ( 10 , -5.0  )  node[above]{{ \\resizebox{ .25\\linewidth}{!}{ \\textcolor{richred}{$\Gamma $} } }}; \n")
    
    if not draw_lset:
        file.write("\\draw[] ( 0.0 , 0.0  )  node[fill=white,above]{{ \\resizebox{ .35\\linewidth}{!}{ \\textcolor{black}{$\mathcal{T}_h $} } }}; \n")
        file.write("\\draw[] ( 0.0 , 7.5  )  node[fill=white,above]{{ \\resizebox{ .35\\linewidth}{!}{ \\textcolor{black}{$\Omega $} } }}; \n")
         
    
    file.write("\\end{tikzpicture} \n") 
    file.write("\\end{document} \n")           
    file.close()

#draw_mesh_simple(name="unfitted-mesh")
#draw_mesh_simple(name="background-mesh",draw_lset=False)


def draw_lset_simple(name="dummy",case_str="geom"):
    ll, ur = (-1.5, -1.5), (1.5, 1.5)
    square = SplineGeometry()
    square.AddRectangle(ll, ur, bc=1)
    mesh = Mesh(square.GenerateMesh(maxh=0.2, quad_dominated=False))
    fes = H1(mesh, order=1, dirichlet=[])
    
    ddx = 10
    rainbow = ["cyan","white"]
    file = open("../plots/meshes/{0}.tex".format(name),"w+")
    file.write("\\documentclass{standalone} \n")
    file.write("\\usepackage{xr} \n")
    file.write("\\usepackage{tikz} \n")
    file.write("\\usepackage{pgfplots} \n")
    file.write("\\usepackage{xcolor} \n")
    file.write("\\usepackage{} \n")   
    file.write("\\usetikzlibrary{shapes,arrows,shadows,snakes,calendar,matrix,spy,backgrounds,folding,calc,positioning,patterns,hobby} \n")
    file.write("\\selectcolormodel{cmyk}  \n")  
    file.write("\\definecolor{orange}{cmyk}{0,0.45,0.91,0} \n")  
    file.write("\\definecolor{brightblue}{cmyk}{0.92,0,0.15,0.05} \n")  
    file.write("\\definecolor{richred}{cmyk}{0,1,0.62,0} \n")  
    file.write("\\definecolor{yellow}{cmyk}{0,0.25,0.95,0} \n")  
    file.write("\\definecolor{navyblue}{cmyk}{1,0.57,0,0.40} \n")  
    file.write("\\begin{document} \n")
    file.write("\\begin{tikzpicture}[scale = 1.0,use Hobby shortcut] \n")  
    file.write("\\pgfresetboundingbox \n")
    file.write("\\path[use as bounding box,draw,black,ultra thick,fill=navyblue, fill opacity=0.8] (-{0},-{0}) rectangle ({0},{0}); \n".format(ddx*1.5) )

    file.write("\\path  \n")
    file.write("(0.0,0.0) coordinate (z0) \n")
    file.write("(2.0,-3.0) coordinate (z1) \n")
    file.write("(6.0,-6.0) coordinate (z2)\n")
    file.write("(10,-2) coordinate (z3)\n")
    file.write("(0.0,10) coordinate (z4)\n")
    file.write("(-10,-2) coordinate (z5)\n")
    file.write("(-6,-6) coordinate (z6)\n")
    file.write("(-2,-3.0) coordinate (z7);\n")
    
    file.write("\draw[richred, line width=2.5mm,closed,fill=white,fill opacity=1.0] (z0) .. (z1) .. (z3) .. (z4) .. (z5) .. (z6) .. (z7); \n")
    file.write("\draw[richred, line width=2.5mm,closed,fill=orange,fill opacity=0.8] (z0) .. (z1) .. (z3) .. (z4) .. (z5) .. (z6) .. (z7); \n")
    
    for el in fes.Elements():
        x_coords = [] 
        y_coords = []
        coords = []
        for vert in el.vertices:
            vx,vy,vz = mesh.ngmesh.Points()[vert.nr+1].p
            coords.append((ddx*vx,ddx*vy))
            x_coords.append(vx)
            y_coords.append(vy)
        x_coords = np.array(x_coords)
        y_coords = np.array(y_coords)
        if len(coords) == 3:
            file.write("\\draw[line width=0.01mm,draw =black]  {1} -- {2} -- {3} -- cycle; \n".format("navyblue",coords[0],coords[1],coords[2] ))
        if len(coords) == 4:
            file.write("\\draw[line width=0.01mm,draw =black] {1} -- {2} -- {3} -- {4} -- cycle; \n".format("white",coords[0],coords[1],coords[2],coords[3]))
    

    file.write("\\draw[] ( 0.0 , 2.5  )  node[fill=white,above]{{   \\resizebox{ .4\\linewidth}{!}{ \\textcolor{black}{$\phi < 0$}}     }}; \n")
    file.write("\\draw[] ( -2 , -9.8  )  node[fill=white,above]{{ \\resizebox{ .4\\linewidth}{!}{ \\textcolor{black}{$\phi >0$ }}      }}; \n")
    file.write("\\draw[] ( 8 , -6.0  )  node[fill=white,above]{{ \\resizebox{ .4\\linewidth}{!}{ \\textcolor{richred}{$\phi = 0 $} } }}; \n")
    
    file.write("\\end{tikzpicture} \n") 
    file.write("\\end{document} \n")           
    file.close()


def fct_on_edge(cf,mesh,va,vb):
    def f(lam):
        return cf(mesh((1-lam)*va[0]+lam*vb[0],(1-lam)*va[1]+lam*vb[1]))
    return f

def FindCut(cf,mesh,va,vb):
    f = fct_on_edge(cf,mesh,va,vb)
    cf_val_va = f(lam=0)
    cf_val_vb = f(lam=1)
    if cf_val_va * cf_val_vb > 1e-14:
        return None
    else:
        sol = optimize.root_scalar(f, bracket=[0, 1], method='brentq')
        print("sol.root = {0}, sol.iterations = {1}, f(root) = {2}".format(sol.root, sol.iterations,f(sol.root) ))
        lam = sol.root
        return (  ddx*((1-lam)*va[0]+lam*vb[0]), ddx*((1-lam)*va[1]+lam*vb[1]) )

def draw_unfitted_space(name="lsetp1.tex",levelset=levelset,dt=NEG):
    ll, ur = (-2, -2), (2, 2)
    square = SplineGeometry()
    square.AddRectangle(ll, ur, bc=1)
    mesh = Mesh(square.GenerateMesh(maxh=0.5, quad_dominated=False))
    fes = H1(mesh, order=1, dirichlet=[])
    
    lsetadap = NoDeformation(mesh, levelset)
    lsetp1 = lsetadap.lset_p1
    ci = CutInfo(mesh, lsetp1)

    hasneg = ci.GetElementsOfType(HASNEG)
    haspos = ci.GetElementsOfType(HASPOS)
    hasif = ci.GetElementsOfType(IF)
    
    print("hasneg =", hasneg)

    rainbow = ["cyan","white"]
    file = open("../plots/meshes/{0}.tex".format(name),"w+")
    file.write("\\documentclass{standalone} \n")
    file.write("\\usepackage{xr} \n")
    file.write("\\usepackage{tikz} \n")
    file.write("\\usepackage{pgfplots} \n")
    file.write("\\usepackage{xcolor} \n")
    file.write("\\usepackage{} \n")   
    file.write("\\usetikzlibrary{shapes,arrows,shadows,snakes,calendar,matrix,spy,backgrounds,folding,calc,positioning,patterns,hobby} \n")
    file.write("\\selectcolormodel{cmyk}  \n")  
    file.write("\\definecolor{orange}{cmyk}{0,0.45,0.91,0} \n")  
    file.write("\\definecolor{brightblue}{cmyk}{0.92,0,0.15,0.05} \n")  
    file.write("\\definecolor{richred}{cmyk}{0,1,0.62,0} \n")  
    file.write("\\definecolor{yellow}{cmyk}{0,0.25,0.95,0} \n")  
    file.write("\\definecolor{navyblue}{cmyk}{1,0.57,0,0.40} \n")  
    file.write("\\begin{document} \n")
    file.write("\\begin{tikzpicture}[scale = 1.0,use Hobby shortcut, spy using outlines= {rectangle, magnification=4, width=19.5cm, height=8cm, connect spies} ] \n")  
    file.write("\\pgfresetboundingbox \n")
    file.write("\\path[use as bounding box,draw,black,ultra thick,fill=white, fill opacity=1.0] (-{0},-{0}) rectangle ({0},{0}); \n".format(ddx*2) )

    
    for el in fes.Elements():
        x_coords = [] 
        y_coords = []
        coords = []

        vert_list = [ ]
        edges = []
        edge_function = [] 
        for vert in el.vertices:
            vx,vy,vz = mesh.ngmesh.Points()[vert.nr+1].p
            vtuple = (vx,vy)
            vert_list.append(vtuple)
        edges = [ [vert_list[0],vert_list[1]] ,  [vert_list[1],vert_list[2]],  [vert_list[2],vert_list[0]] ]
        cut_on_edge = [ ] 
        for edge in edges: 
            result = FindCut( lsetp1 , mesh, edge[0], edge[1])
            if result: 
                cut_on_edge.append(result)
        print("cut_on_edge = ", cut_on_edge)


        
        for vert in el.vertices:
            vx,vy,vz = mesh.ngmesh.Points()[vert.nr+1].p
            coords.append((ddx*vx,ddx*vy))
            x_coords.append(vx)
            y_coords.append(vy)
        x_coords = np.array(x_coords)
        y_coords = np.array(y_coords)
        if dt == NEG:
            if len(coords) == 3 and hasneg[el.nr] :
                file.write("\\draw[line width=0.01mm,draw =black,fill=orange,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
            else:
                file.write("\\draw[line width=0.01mm,draw =black,fill=white,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
        if dt == POS:
            if len(coords) == 3 and haspos[el.nr] :
                file.write("\\draw[line width=0.01mm,draw =black,fill=navyblue,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
            else:
                file.write("\\draw[line width=0.01mm,draw =black,fill=white,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
        if dt == IF:
            if len(coords) == 3 and hasif[el.nr] :
                file.write("\\draw[line width=0.01mm,draw =black,fill=yellow,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
            else:
                file.write("\\draw[line width=0.01mm,draw =black,fill=white,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))

        if cut_on_edge != []:
            #vx,vy = cut_on_edge[0][0],cut_on_edge[0][1]
            #wx,wy = cut_on_edge[1][0],cut_on_edge[1][1]
            if dt == IF:
                file.write("\\draw[line width=2mm,draw=green!70!black]  {0} -- {1}; \n".format(cut_on_edge[0] , cut_on_edge[1] ))
            if dt == NEG:
                file.write("\\draw[line width=1mm,draw=green!70!black]  {0} -- {1}; \n".format(cut_on_edge[0] , cut_on_edge[1] ))

    #file.write("\\draw[] ( 0.0 , 2.5  )  node[fill=white,above]{{   \\resizebox{ .4\\linewidth}{!}{ \\textcolor{black}{$\phi < 0$}}     }}; \n")
    #file.write("\\draw[] ( -2 , -9.8  )  node[fill=white,above]{{ \\resizebox{ .4\\linewidth}{!}{ \\textcolor{black}{$\phi >0$ }}      }}; \n")
    #file.write("\\draw[] ( 0.0 , -10.0  )  node[above]{{ \\resizebox{ .25\\linewidth}{!}{ \\textcolor{richred}{$\Gamma $} } }}; \n")
    if dt==NEG or dt==POS: 
        file.write("\\begin{{axis}}[view={{0}}{{90}},  anchor=origin,  disabledatascaling, at={{(0pt,0pt)}}, x={0}cm,y={0}cm,z={0}cm, hide axis ] \n".format(ddx))
        #file.write("\\addplot3 [ultra thick, \n")
        if dt == NEG:
            file.write("\\addplot3[line width=1mm,  \n") 
        else:
            file.write("\\addplot3[line width=2mm,  \n") 
        #file.write(" contour lua={levels={0},labels=false, draw color=richred}, samples=200 ] {x^4 + y^4 - 1}; \n")
        file.write(" contour lua={levels={0},labels=false, draw color=richred}, samples=200 ] {x^2 + y^2 - 1}; \n")
        file.write(" \end{axis} \n ")
    if dt == NEG:
        file.write("\\draw[] ( 0.0 , 0.0  )  node[fill=white,above]{{ \\resizebox{ .35\\linewidth}{!}{ \\textcolor{orange}{$\mathcal{T}_h^1 $} } }}; \n")
        file.write("\\draw[] ( -3.5 , -10.45  )  node[above]{{ \\resizebox{ .12\\linewidth}{!}{ \\textcolor{richred}{$ \Gamma$} } }}; \n")
        file.write("\\draw[] ( -0.7 , -9.6  )  node[above]{{ \\resizebox{ .2\\linewidth}{!}{ \\textcolor{green!70!black}{ $\Gamma^{\mathrm{lin}} $  } } }}; \n")
        file.write("\\coordinate (spypoint) at (-1.75,-9.5);  \n" )
        file.write("\\coordinate (magnifyglass) at (10,-15.75); \n")
        file.write("\\spy [black] on (spypoint) in node[fill=white] at (magnifyglass); \n")
    if dt == POS:
        file.write("\\draw[] ( -12 , -12  )  node[fill=white,above]{{ \\resizebox{ .35\\linewidth}{!}{ \\textcolor{navyblue}{$\mathcal{T}_h^2 $} } }}; \n")
        file.write("\\draw[] ( 0.0 , -13.25  )  node[above]{{ \\resizebox{ .35\\linewidth}{!}{ \\textcolor{richred}{$\Gamma $} } }}; \n")
    if dt == IF:
        file.write("\\draw[] ( -9.75 , 0  )  node[fill=white,above]{{ \\resizebox{ .35\\linewidth}{!}{ \\textcolor{yellow}{$\mathcal{T}_h^{\Gamma} $} } }}; \n")
        file.write("\\draw[] ( 0 , 13.5  )  node[fill=white,above]{{ \\resizebox{ .35\\linewidth}{!}{ \\textcolor{black}{$\mathcal{T}_h $} } }}; \n")
        file.write("\\draw[] ( -0.7 , -9.25  )  node[above]{{ \\resizebox{ .55\\linewidth}{!}{ \\textcolor{green!70!black}{ $\Gamma^{\mathrm{lin}} $  } } }}; \n")

    file.write("\\end{tikzpicture} \n") 
    file.write("\\end{document} \n")           
    file.close()

#draw_unfitted_space(name="unfitted-neg",levelset=levelset,dt=NEG)
#draw_unfitted_space(name="unfitted-if",levelset=levelset,dt=IF)
#draw_unfitted_space(name="unfitted-pos",levelset=levelset,dt=POS)


def draw_GeomStab(name="GeomStab",levelset=levelset,dt=NEG):
    ll, ur = (-2, -2), (2, 2)
    square = SplineGeometry()
    square.AddRectangle(ll, ur, bc=1)
    mesh = Mesh(square.GenerateMesh(maxh=0.25, quad_dominated=False))
    fes = H1(mesh, order=1, dirichlet=[])
    
    lsetadap = NoDeformation(mesh, levelset)
    lsetp1 = lsetadap.lset_p1
    ci = CutInfo(mesh, lsetp1)

    hasneg = ci.GetElementsOfType(HASNEG)
    haspos = ci.GetElementsOfType(HASPOS)
    hasif = ci.GetElementsOfType(IF)
    
    print("hasneg =", hasneg)


    Vaux = L2(mesh, order=0, dirichlet=[],dgjumps=True)
    support_dfm  = GridFunction(Vaux)
    lsetadap_aux = LevelSetMeshAdaptation(mesh, order=2, levelset=levelset)
    support_dfm.Set(IfPos(lsetadap_aux.deform.Norm(),1.0,0.0))
    for i in range(Vaux.ndof):
        if( support_dfm.vec[i]) > 1e-14:
            support_dfm.vec[i] = 1.0

    rainbow = ["cyan","white"]
    file = open("../plots/meshes/{0}.tex".format(name),"w+")
    file.write("\\documentclass{standalone} \n")
    file.write("\\usepackage{xr} \n")
    file.write("\\usepackage{tikz} \n")
    file.write("\\usepackage{pgfplots} \n")
    file.write("\\usepackage{xcolor} \n")
    file.write("\\usepackage{} \n")   
    file.write("\\usetikzlibrary{shapes,arrows,shadows,snakes,calendar,matrix,spy,backgrounds,folding,calc,positioning,patterns,hobby} \n")
    file.write("\\selectcolormodel{cmyk}  \n")  
    file.write("\\definecolor{orange}{cmyk}{0,0.45,0.91,0} \n")  
    file.write("\\definecolor{brightblue}{cmyk}{0.92,0,0.15,0.05} \n")  
    file.write("\\definecolor{richred}{cmyk}{0,1,0.62,0} \n")  
    file.write("\\definecolor{yellow}{cmyk}{0,0.25,0.95,0} \n")  
    file.write("\\definecolor{navyblue}{cmyk}{1,0.57,0,0.40} \n")  
    file.write("\\begin{document} \n")
    file.write("\\begin{tikzpicture}[scale = 1.0,use Hobby shortcut, spy using outlines= {rectangle, magnification=4, width=19.5cm, height=8cm, connect spies} ] \n")  
    file.write("\\pgfresetboundingbox \n")
    file.write("\\path[use as bounding box,draw,black,ultra thick,fill=white, fill opacity=1.0] (-{0},-{0}) rectangle ({0},{0}); \n".format(ddx*2) )

    
    for el in fes.Elements():
        x_coords = [] 
        y_coords = []
        coords = []

        vert_list = [ ]
        edges = []
        edge_function = [] 
        for vert in el.vertices:
            vx,vy,vz = mesh.ngmesh.Points()[vert.nr+1].p
            vtuple = (vx,vy)
            vert_list.append(vtuple)
        edges = [ [vert_list[0],vert_list[1]] ,  [vert_list[1],vert_list[2]],  [vert_list[2],vert_list[0]] ]
        cut_on_edge = [ ] 
        for edge in edges: 
            result = FindCut( lsetp1 , mesh, edge[0], edge[1])
            if result: 
                cut_on_edge.append(result)
        print("cut_on_edge = ", cut_on_edge)


        cent_x = 0.0 
        cent_y = 0.0
        for vert in el.vertices:
            vx,vy,vz = mesh.ngmesh.Points()[vert.nr+1].p
            coords.append((ddx*vx,ddx*vy))
            x_coords.append(vx)
            y_coords.append(vy)
            cent_x += vx 
            cent_y += vy
        x_coords = np.array(x_coords)
        y_coords = np.array(y_coords)
        cent_x *= 1/3
        cent_y *= 1/3
        if dt == NEG:
            if len(coords) == 3 and hasneg[el.nr] :
                file.write("\\draw[line width=0.01mm,draw =black,fill=orange,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
            else:
                file.write("\\draw[line width=0.01mm,draw =black,fill=white,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
        if dt == POS:
            if len(coords) == 3 and haspos[el.nr] :
                file.write("\\draw[line width=0.01mm,draw =black,fill=navyblue,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
            else:
                file.write("\\draw[line width=0.01mm,draw =black,fill=white,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
        if dt == IF:
            if support_dfm(mesh(cent_x,cent_y)) > 1e-14: 
                file.write("\\draw[line width=0.01mm,draw =black,fill=cyan,fill opacity=0.3]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
            else:
                file.write("\\draw[line width=0.01mm,draw =black,fill=white,fill opacity=1.0]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))

            #if len(coords) == 3 and hasif[el.nr] :
            #    file.write("\\draw[line width=0.01mm,draw =black,fill=yellow,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
            #else:
            #    file.write("\\draw[line width=0.01mm,draw =black,fill=white,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))

        if cut_on_edge != []:
            #vx,vy = cut_on_edge[0][0],cut_on_edge[0][1]
            #wx,wy = cut_on_edge[1][0],cut_on_edge[1][1]
            if dt == IF:
                file.write("\\draw[line width=2mm,draw=green!70!black]  {0} -- {1}; \n".format(cut_on_edge[0] , cut_on_edge[1] ))
            if dt == NEG:
                file.write("\\draw[line width=1mm,draw=green!70!black]  {0} -- {1}; \n".format(cut_on_edge[0] , cut_on_edge[1] ))

    #file.write("\\draw[] ( 0.0 , 2.5  )  node[fill=white,above]{{   \\resizebox{ .4\\linewidth}{!}{ \\textcolor{black}{$\phi < 0$}}     }}; \n")
    #file.write("\\draw[] ( -2 , -9.8  )  node[fill=white,above]{{ \\resizebox{ .4\\linewidth}{!}{ \\textcolor{black}{$\phi >0$ }}      }}; \n")
    #file.write("\\draw[] ( 0.0 , -10.0  )  node[above]{{ \\resizebox{ .25\\linewidth}{!}{ \\textcolor{richred}{$\Gamma $} } }}; \n")
    if dt==NEG or dt==POS: 
        file.write("\\begin{{axis}}[view={{0}}{{90}},  anchor=origin,  disabledatascaling, at={{(0pt,0pt)}}, x={0}cm,y={0}cm,z={0}cm, hide axis ] \n".format(ddx))
        #file.write("\\addplot3 [ultra thick, \n")
        if dt == NEG:
            file.write("\\addplot3[line width=1mm,  \n") 
        else:
            file.write("\\addplot3[line width=2mm,  \n") 
        #file.write(" contour lua={levels={0},labels=false, draw color=richred}, samples=200 ] {x^4 + y^4 - 1}; \n")
        file.write(" contour lua={levels={0},labels=false, draw color=richred}, samples=200 ] {x^2 + y^2 - 1}; \n")
        file.write(" \end{axis} \n ")
    if dt == NEG:
        file.write("\\draw[] ( 0.0 , 0.0  )  node[fill=white,above]{{ \\resizebox{ .35\\linewidth}{!}{ \\textcolor{orange}{$\mathcal{T}_h^1 $} } }}; \n")
        file.write("\\draw[] ( -3.5 , -10.45  )  node[above]{{ \\resizebox{ .12\\linewidth}{!}{ \\textcolor{richred}{$ \Gamma$} } }}; \n")
        file.write("\\draw[] ( -0.7 , -9.6  )  node[above]{{ \\resizebox{ .2\\linewidth}{!}{ \\textcolor{green!70!black}{ $\Gamma^{\mathrm{lin}} $  } } }}; \n")
        file.write("\\coordinate (spypoint) at (-1.75,-9.5);  \n" )
        file.write("\\coordinate (magnifyglass) at (10,-15.75); \n")
        file.write("\\spy [black] on (spypoint) in node[fill=white] at (magnifyglass); \n")
    if dt == POS:
        file.write("\\draw[] ( -12 , -12  )  node[fill=white,above]{{ \\resizebox{ .35\\linewidth}{!}{ \\textcolor{navyblue}{$\mathcal{T}_h^2 $} } }}; \n")
        file.write("\\draw[] ( 0.0 , -13.25  )  node[above]{{ \\resizebox{ .35\\linewidth}{!}{ \\textcolor{richred}{$\Gamma $} } }}; \n")
    if dt == IF:
        #file.write("\\draw[] ( -9.75 , 0  )  node[fill=white,above]{{ \\resizebox{ .35\\linewidth}{!}{ \\textcolor{yellow}{$\mathcal{T}_h^{\Gamma} $} } }}; \n")
        #file.write("\\draw[] ( 0 , 13.5  )  node[fill=white,above]{{ \\resizebox{ .35\\linewidth}{!}{ \\textcolor{black}{$\mathcal{T}_h $} } }}; \n")
        file.write("\\draw[fill=white] ( -0.2 , -8.25  )  node[fill=white,above]{{ \\resizebox{ .65\\linewidth}{!}{ \\textcolor{green!70!black}{ $\Gamma^{\mathrm{lin}} $  } } }}; \n")
        file.write("\\draw[fill=white] ( -0.2 , 11.75  )  node[fill=white,above]{{ \\resizebox{ 1.25\\linewidth}{!}{ \\textcolor{cyan}{ $  \mathcal{T}^{\Gamma,1}_{+} \cup \mathcal{T}^{\Gamma,2}_{+} $  } } }}; \n")

    file.write("\\end{tikzpicture} \n") 
    file.write("\\end{document} \n")           
    file.close()

#draw_GeomStab(name="GeomStab",levelset=levelset,dt=IF)


def draw_Stab(name="Stab-sketch.tex",levelset=levelset,dt=NEG):
    ll, ur = (-2, -2), (2, 2)
    square = SplineGeometry()
    square.AddRectangle(ll, ur, bc=1)
    mesh = Mesh(square.GenerateMesh(maxh=0.465, quad_dominated=False))
    fes = H1(mesh, order=1, dirichlet=[])
    
    lsetadap = NoDeformation(mesh, levelset)
    lsetp1 = lsetadap.lset_p1
    ci = CutInfo(mesh, lsetp1)

    hasneg = ci.GetElementsOfType(HASNEG)
    haspos = ci.GetElementsOfType(HASPOS)
    hasif = ci.GetElementsOfType(IF)
    
    print("hasneg =", hasneg)

    rainbow = ["cyan","white"]
    file = open("{0}.tex".format(name),"w+")
    file.write("\\documentclass{standalone} \n")
    file.write("\\usepackage{xr} \n")
    file.write("\\usepackage{tikz} \n")
    file.write("\\usepackage{pgfplots} \n")
    file.write("\\usepackage{xcolor} \n")
    file.write("\\usepackage{} \n")   
    file.write("\\usetikzlibrary{shapes,arrows,shadows,snakes,calendar,matrix,spy,backgrounds,folding,calc,positioning,patterns,hobby} \n")
    file.write("\\selectcolormodel{cmyk}  \n")  
    file.write("\\definecolor{orange}{cmyk}{0,0.45,0.91,0} \n")  
    file.write("\\definecolor{brightblue}{cmyk}{0.92,0,0.15,0.05} \n")  
    file.write("\\definecolor{richred}{cmyk}{0,1,0.62,0} \n")  
    file.write("\\definecolor{yellow}{cmyk}{0,0.25,0.95,0} \n")  
    file.write("\\definecolor{navyblue}{cmyk}{1,0.57,0,0.40} \n")  
    file.write("\\begin{document} \n")
    file.write("\\begin{tikzpicture}[scale = 1.0,use Hobby shortcut, spy using outlines= {rectangle, magnification=4, width=19.5cm, height=8cm, connect spies} ] \n")  
    file.write("\\pgfresetboundingbox \n")
    file.write("\\path[use as bounding box,draw,black,ultra thick,fill=white, fill opacity=1.0] (-{0},-{0}) rectangle ({0},{0}); \n".format(ddx*2) )

    
    for el in fes.Elements():
        x_coords = [] 
        y_coords = []
        coords = []

        vert_list = [ ]
        edges = []
        edge_function = [] 
        for vert in el.vertices:
            vx,vy,vz = mesh.ngmesh.Points()[vert.nr+1].p
            vtuple = (vx,vy)
            vert_list.append(vtuple)
        edges = [ [vert_list[0],vert_list[1]] ,  [vert_list[1],vert_list[2]],  [vert_list[2],vert_list[0]] ]
        cut_on_edge = [ ] 
        for edge in edges: 
            result = FindCut( lsetp1 , mesh, edge[0], edge[1])
            if result: 
                cut_on_edge.append(result)
        print("cut_on_edge = ", cut_on_edge)
        
        centroid_x = 0
        centroid_y = 0
        for vert in el.vertices:
            vx,vy,vz = mesh.ngmesh.Points()[vert.nr+1].p
            coords.append((ddx*vx,ddx*vy))
            centroid_x += ddx*vx
            centroid_y += ddx*vy
            x_coords.append(vx)
            y_coords.append(vy)
        x_coords = np.array(x_coords)
        y_coords = np.array(y_coords)
        centroid_x *= 1/3
        centroid_y *= 1/3
        if centroid_x > -1.5*ddx and centroid_x <  0.1*ddx and centroid_y > -0.6*ddx and  centroid_y < 0.9*ddx:
            file.write("\\draw[line width=0.01mm,draw =black,fill=white,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
        #if dt == NEG:
        #    if len(coords) == 3 and hasneg[el.nr] :
        #        file.write("\\draw[line width=0.01mm,draw =black,fill=orange,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
        #    else:
        #        file.write("\\draw[line width=0.01mm,draw =black,fill=white,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
        #if dt == POS:
        #    if len(coords) == 3 and haspos[el.nr] :
        #        file.write("\\draw[line width=0.01mm,draw =black,fill=navyblue,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
        #    else:
        #        file.write("\\draw[line width=0.01mm,draw =black,fill=white,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
        #if dt == IF:
        #    if len(coords) == 3 and hasif[el.nr] :
        #        file.write("\\draw[line width=0.01mm,draw =black,fill=yellow,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))
        #    else:
        #        file.write("\\draw[line width=0.01mm,draw =black,fill=white,fill opacity=0.8]  {0} -- {1} -- {2} -- cycle; \n".format(coords[0],coords[1],coords[2] ))

        if cut_on_edge != []:
            if centroid_x > -1.45*ddx and centroid_x <  0.1*ddx and centroid_y > -0.95*ddx and  centroid_y < 1.0*ddx:
                file.write("\\draw[line width=2mm,draw=green!70!black]  {0} -- {1}; \n".format(cut_on_edge[0] , cut_on_edge[1] ))
            #vx,vy = cut_on_edge[0][0],cut_on_edge[0][1]
            #wx,wy = cut_on_edge[1][0],cut_on_edge[1][1]

    file.write("\\end{tikzpicture} \n") 
    file.write("\\end{document} \n")           
    file.close()

draw_Stab(name="Stab-sketch",levelset=levelset,dt=NEG)
