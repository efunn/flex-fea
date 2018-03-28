import dolfin as df
from mshr import *
from math import sin, cos, pi
import numpy as np

class ParallelFlex(object):
    def __init__(self):
        # geometry constants
        self.end_thk = 5*1e-3
        self.wall_gap = 2.4*1e-3
        self.mirror_width = 5*1e-3
        self.press_len = 2*1e-3
        self.press_width = 0.5*1e-3
        self.total_thk = 7.5e-3
        self.top_key_height = 3e-3
        self.support_poly_gap = 1e-3

        # PDE params (need to translate to FEniCS from Matlab)
        self.pde_e = 1300e6 #2310e6 # Pa (reported between 1000 and 2000 MPa) # 1300 MPa from other datasheet?
        self.pde_gnu = .408 # Poisson's ratio of the material: 0.408
        self.pde_f = np.array([0, 0]).T # No body forces
        G = self.pde_e/(2.*(1+self.pde_gnu))
        # mu = 2*G*self.pde_gnu/(1-self.pde_gnu) # plane stress
        mu = 2*G*self.pde_gnu/(1-2*self.pde_gnu) # plane strain
        self.pde_c = np.array([2*G+mu, 0, G,   0, G, mu, 0,  G, 0, 2*G+mu]).T
        # TODO: create (vector) function space
        # V = VectorFunctionSpace(self.mesh, "CG", 1)
        # create boundary conditions:
        # class Border(SubDomain):
        #     def inside(self, x, on_boundary):
        #         return x[1]==0.0
        # bc = DirichletBC(V, 0., border) 

    def set_geometry(self, wall_thk=1, wall_len=15,
            theta=45, pinch_len=40, press_pos=0, resolution=45):
        len_ratio = 1.0;
        thk_ratio = 1.0;
        taper_ratio = 1.0;
        self.wall_thk = wall_thk*1e-3;
        self.wall_len = wall_len*1e-3;
        self.pinch_len = pinch_len*1e-3;
        self.press_pos = press_pos*1e-3;
        self.angle_wall_len = self.wall_len*len_ratio;
        self.theta = theta;

        # dependents
        self.angle_len = self.angle_wall_len*sin(deg2rad(theta));
        self.total_width = 4*self.wall_thk+2*self.wall_gap+2*self.angle_len+self.mirror_width;
        self.total_len = 2*self.end_thk+2*self.wall_len+self.pinch_len;

        # csg geometry for basic geometry
        body = Rectangle(df.Point(0,0),df.Point(self.total_len,self.total_width))

        x_off = self.wall_len+self.pinch_len;
        y_off = 2*self.wall_thk+self.wall_gap+2*self.angle_len+self.mirror_width;
        cr_p0 = [self.end_thk,self.wall_thk]
        cr_p1 = [self.end_thk+self.wall_len,self.wall_thk+self.wall_gap]
        corner_rect_1 = Rectangle(df.Point(cr_p0[0],cr_p0[1]),df.Point(cr_p1[0],cr_p1[1]))
        corner_rect_2 = Rectangle(df.Point(cr_p0[0]+x_off,cr_p0[1]),df.Point(cr_p1[0]+x_off,cr_p1[1]))
        corner_rect_3 = Rectangle(df.Point(cr_p0[0],cr_p0[1]+y_off),df.Point(cr_p1[0],cr_p1[1]+y_off))
        corner_rect_4 = Rectangle(df.Point(cr_p0[0]+x_off,cr_p0[1]+y_off),df.Point(cr_p1[0]+x_off,cr_p1[1]+y_off))

        x1 = self.end_thk;
        x2 = x1+self.wall_len;
        x3 = x2+self.angle_wall_len*cos(deg2rad(self.theta));
        x4 = self.end_thk+self.wall_len+self.pinch_len;
        x5 = x4+self.angle_wall_len*cos(deg2rad(self.theta));
        x6 = x4+self.wall_len;
        x7 = x2+self.wall_thk/sin(deg2rad(self.theta));#*taper_ratio;
        x8 = x3+self.wall_thk/sin(deg2rad(self.theta));
        x9 = x5-self.wall_thk/sin(deg2rad(self.theta))*thk_ratio;
        x10 = x4-self.wall_thk/sin(deg2rad(self.theta))*thk_ratio;#*taper_ratio;
        y1 = 2*self.wall_thk+self.wall_gap;
        y2 = y1+2*self.angle_len+self.mirror_width;
        y3 = y2-self.angle_len;
        y4 = y3-self.mirror_width;
        left_poly = Polygon([df.Point(x1,y1),df.Point(x1,y2),df.Point(x2,y2),
            df.Point(x3,y3),df.Point(x3,y4),df.Point(x2,y1)][::-1])
        right_poly = Polygon([df.Point(x4,y1),df.Point(x5,y4),df.Point(x5,y3),
            df.Point(x4,y2),df.Point(x6,y2),df.Point(x6,y1)][::-1])
        mid_poly_1 = Polygon([df.Point(x7,y1),df.Point(x8,y4),df.Point(x9,y4),df.Point(x10,y1)][::-1])
        poly_y_off = self.angle_len+self.mirror_width;
        mid_poly_2 = Polygon([df.Point(x7,y4+poly_y_off),df.Point(x8,y1+poly_y_off),
            df.Point(x9,y1+poly_y_off),df.Point(x10,y4+poly_y_off)])

        press_x1 = .5*self.total_len-.5*self.pinch_len;
        press_x2 = press_x1+self.pinch_len;
        press_x3 = .5*self.total_len-.5*self.press_len+self.press_pos;
        press_x4 = press_x3+self.press_len;
        press_y1 = 0;
        press_y2 = -self.press_width;
        top_key_y1 = self.total_width;
        top_key_y2 = self.total_width+self.top_key_height;
        press_y3 = top_key_y2;
        press_y4 = press_y3+self.press_width;
        bottom_press_rect = Rectangle(df.Point(press_x1,press_y1),df.Point(press_x2,press_y2))
        top_press_rect = Rectangle(df.Point(press_x3,press_y3),df.Point(press_x4,press_y4))
        # press_support poly
        sup_x2 = x7+self.support_poly_gap/sin(deg2rad(self.theta));
        sup_x1 = sup_x2+(self.angle_wall_len-self.support_poly_gap/sin(deg2rad(self.theta)))*cos(deg2rad(self.theta));
        sup_x3 = x10-self.support_poly_gap/sin(deg2rad(self.theta));
        sup_x4 = sup_x3-3e-3;
        sup_y1 = y1+poly_y_off+self.support_poly_gap;
        sup_y2 = y2;
        press_support_poly = Polygon([df.Point(sup_x1,sup_y1),df.Point(sup_x2,sup_y2),
            df.Point(sup_x3,sup_y2),df.Point(sup_x4,y1)][::-1])
        top_key_rect = Rectangle(df.Point(sup_x2,top_key_y1),df.Point(sup_x3,top_key_y2))

        self.domain = (body+bottom_press_rect+top_press_rect+top_key_rect
            -corner_rect_1-corner_rect_2-corner_rect_3-corner_rect_4 
            -left_poly-right_poly-mid_poly_1-mid_poly_2)#+press_support_poly)
        self.mesh = generate_mesh(self.domain,resolution)

    def show_geometry(self):
        df.plot(self.mesh,interactive=True)

def deg2rad(degrees):
    return degrees*pi/180.

if __name__ == '__main__':
    flex = ParallelFlex()
    flex.set_geometry()
    flex.show_geometry()
