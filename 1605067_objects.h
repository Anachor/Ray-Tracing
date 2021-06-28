#ifndef __1605067_OBJECTS_H__
#define __1605067_OBJECTS_H__

#include "1605067_geometry.h"
#include<tuple>

namespace Objects {
    const double INF = 1e100;
    using Geometry::Point;
    struct Color{
        double r, g, b; //Between [0, 1]

        friend std::istream &operator >> (std::istream &is, Color &o) { 
            return is >> o.r >> o.g >> o.b; 
        }
        friend std::ostream &operator << (std::ostream &os, const Color &o) { 
            return os << o.r <<" "<< o.g <<" "<< o.b; 
        }
    };

    struct Ray{
        Point a, dir;
        Ray(Point a, Point dir) : a(a), dir(Geometry::unit(dir)) {}
    };


    struct Reflection_coefficients {
        double ambient, diffuse, specular, recursive, exponent;
    };


    struct Light {
        
        Point position;
        Color color;

        void draw() {
            glPointSize( 6.0 );
            glColor3f(color.r, color.g, color.b);
            glBegin(GL_POINTS);{
                glVertex3f(position.x, position.y, position.z);

            }glEnd();
        }
    };

    struct Object{
        Color color;
        Reflection_coefficients coefficients;

        Object( Color color, Reflection_coefficients coefficients)
            : color(color), coefficients(coefficients)
        {}

        Object() {}

        virtual void draw() = 0;
        //virtual Geometry::Point normal() {}
        virtual std::pair<double, Color> intersect(Ray r) = 0;
    };

    struct Floor : Object {
        double floorWidth;
        double tileWidth;
        Color darkcolor, lightcolor;

        Floor(double floorWidth, double tileWidth, Color darkcolor, Color lightcolor) 
            : floorWidth(floorWidth), tileWidth(tileWidth), darkcolor(darkcolor), lightcolor(lightcolor)
        {}
        std::pair<double, Color> intersect(Ray r) override { 
            double len = Geometry::length(r.dir);
            double d = -r.a.z/r.dir.z;
            Point P = r.a + r.dir*d;
            if (!(abs(P.x) <= floorWidth/2 && abs(P.y) <= floorWidth/2))    return {INF, {}};
            int i = (P.x+floorWidth/2)/tileWidth;
            int j = (P.y+floorWidth/2)/tileWidth;
            if (i%2 == j%2)     return {d, darkcolor};
            else                return {d, lightcolor};
        };
        void draw(){
            assert(floorWidth > 0 && tileWidth > 0 && floorWidth/tileWidth <= 1000);
            for (auto [x, i] = std::make_tuple(-floorWidth/2, 0); x<=floorWidth/2; x += tileWidth, i++) {
                for (auto [y, j] = std::make_tuple(-floorWidth/2, 0); y<=floorWidth/2; y += tileWidth, j++) {
                    Color c = ( i%2 == j%2 ? darkcolor : lightcolor);
                    glColor3f(c.r, c.g, c.b);


                    glBegin(GL_QUADS);{
                        glVertex3f(x, y, 0);
                        glVertex3f(x+tileWidth,  y, 0);
                        glVertex3f(x+tileWidth, y+tileWidth, 0);
                        glVertex3f(x, y+tileWidth, 0);

                    }glEnd();
                }   
            }
        }
    };


    struct Sphere : Object{
        Point center;
        double radius;

        const static int stacks = 100, slices = 100;

        Sphere(Point center, double radius) : center(center), radius(radius) {
        }

        std::pair<double, Color> intersect(Ray r) override {
            // std::cout<<"Ray: "<<r.a<<" --> "<<r.dir<<"  ----  Sphere: O = "<<center<<", r = "<<radius<<std::endl;

            Geometry::Line L = {r.a, r.dir};
            double OX = Geometry::distancePointLine(center, L);
            double OA = length(r.a-center); 
            
            double D = sqrt(OA*OA-OX*OX);
            double d = sqrt(radius*radius-OX*OX);

            Point Z1 = r.a + unit(r.dir)*D;
            Point Z2 = r.a - unit(r.dir)*D;
            if (Geometry::length(Z2 - center) < Geometry::length(Z1 - center))    D = -D;

            
            if (OX > radius)  return {INF, {}};
            else if (D+d < 0)            return {INF, {}};
            else if (D-d >= 0)      return {D-d, color};
            else                    return {D+d, color};
        }
        
        void draw() override {
            // write codes for drawing sphere
            using Geometry::Point, Geometry::PI;
            Point points[stacks+1][slices+1];

            //generate points
            for(int i=0;i<=stacks;i++) {
                double h=radius*sin(((double)i/(double)stacks)*(PI/2));
                double r=radius*cos(((double)i/(double)stacks)*(PI/2));
                for(int j=0;j<=slices;j++) {
                    points[i][j].x=r*cos(((double)j/(double)slices)*2*PI);
                    points[i][j].y=r*sin(((double)j/(double)slices)*2*PI);
                    points[i][j].z=h;
                }
            }
            //draw quads using generated points
            for(int i=0;i<stacks;i++) {
                for(int j=0;j<slices;j++) {
                    glBegin(GL_QUADS);{
                        glColor3f(color.r, color.g, color.b);

                        //upper hemisphere
                        glVertex3f(center.x+points[i][j].x,center.y+points[i][j].y,center.z+points[i][j].z);
                        glVertex3f(center.x+points[i][j+1].x,center.y+points[i][j+1].y,center.z+points[i][j+1].z);
                        glVertex3f(center.x+points[i+1][j+1].x,center.y+points[i+1][j+1].y,center.z+points[i+1][j+1].z);
                        glVertex3f(center.x+points[i+1][j].x,center.y+points[i+1][j].y,center.z+points[i+1][j].z);
                        //lower hemisphere
                        glVertex3f(center.x+points[i][j].x,center.y+points[i][j].y,center.z-points[i][j].z);
                        glVertex3f(center.x+points[i][j+1].x,center.y+points[i][j+1].y,center.z-points[i][j+1].z);
                        glVertex3f(center.x+points[i+1][j+1].x,center.y+points[i+1][j+1].y,center.z-points[i+1][j+1].z);
                        glVertex3f(center.x+points[i+1][j].x,center.y+points[i+1][j].y,center.z-points[i+1][j].z);
                    }glEnd();
                }
            }
        }
    };


    struct Triangle : Object{
        Point A, B, C;
        
        Triangle(Point A, Point B, Point C) : A(A), B(B), C(C) {
        }
        std::pair<double, Color> intersect(Ray r) override { 
            using namespace Geometry;

            Plane p(A, B, C);
            double d = p.intersectLine({r.a, r.dir})/length(r.dir);
            if (d < 0)      return {INF, {}};

            Point O = r.a + unit(r.dir)*d;
            double ABC = length(cross(A-B, B-C));
            double OAB = length(cross(A-O, B-O));
            double OBC = length(cross(B-O, C-O));
            double OCA = length(cross(C-O, A-O));

            if (Geometry::dcmp((OAB + OBC + OCA)/ABC - 1) <= 0)     return {d, color};
            else                            return {INF, {}};
        };
        void draw() override {
            // write codes for drawing triangle
            glColor3f(color.r, color.g, color.b);
            glBegin(GL_TRIANGLES);{
                glVertex3f(A.x, A.y, A.z);
                glVertex3f(B.x, B.y, B.z);
                glVertex3f(C.x, C.y, C.z);

            }glEnd();
        }
    };

    std::vector<double> solveQuadratic(double A, double B, double C) {
        using Geometry::EPS;
        if (abs(A) < EPS) {
            if (abs(B) < EPS)   return {};
            return {-C/B};
        }
        else {
            double D = B*B-4*A*C;
            if (D < -EPS)   return {};
            else if (abs(D) < EPS)    return {-B/(2*A)};
            else return { (-B-sqrt(D)/(2*A)), (-B+sqrt(D)/(2*A))};
        }
    }

    struct GeneralQuadraticSurface : Object{
        ///Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fzx + Gx + Hy + Iz + J
        double A, B, C, D, E, F, G, H, I, J;
        Geometry::Point reference_point;
        double length, width, height;

        bool inRange(double x, double len, double X) {
            return len == 0 || (x <= X && X <= x+len);
        }

        bool inCube(Point Z) {
            return inRange(reference_point.x, length, Z.x) && 
            inRange(reference_point.y, width, Z.y) &&
            inRange(reference_point.z, height, Z.z);
        }
        
        GeneralQuadraticSurface(double A, double B, double C, double D, double E,
                                double F, double G, double H, double I, double J,
                                Geometry::Point reference_point,
                                double length, double width, double height) 
            : A(A), B(B), C(C), D(D), E(E), F(F), G(G), H(H), I(I), J(J),
              reference_point(reference_point), length(length), width(width), height(height)
        {}
        std::pair<double, Color> intersect(Ray r) override { 
            
            double qA =  A*r.dir.x*r.dir.x + B*r.dir.y*r.dir.y + C*r.dir.z*r.dir.z 
                       + D*r.dir.x*r.dir.y + E*r.dir.x*r.dir.z + F*r.dir.y*r.dir.z;

            double qB =  2*A*r.a.x*r.dir.x + 2*B*r.a.y*r.dir.y + 2*C*r.a.z*r.dir.z 
                        + D*(r.a.x*r.dir.y + r.a.y*r.dir.x)
                        + E*(r.a.z*r.dir.x + r.a.x*r.dir.z)
                        + F*(r.a.y*r.dir.z + r.a.z*r.dir.y)
                        + G*r.dir.x + H*r.dir.y + I*r.dir.z;

            double qC =   A*r.a.x*r.a.x + B*r.a.y*r.a.y + C*r.a.z*r.a.z
                        + D*r.a.x*r.a.y + E*r.a.x*r.a.z + F*r.a.y*r.a.z 
                        + G*r.a.x + H*r.a.y + I*r.a.z + J;
                

            double closest = INF, len = Geometry::length(r.dir);
            for (double d: solveQuadratic(qA, qB, qC)) {
                Point p = r.a + r.dir*d;
                if (!inCube(p)) continue;
                d = d*len;
                if (d < closest && d >= 0) {
                    closest = d;
                }
            }

            return {closest, color};
        };
        void draw() override {
            //std::cout<<"Cannot draw general quadratic surface"<<std::endl;
        }
    };

}

#endif // __1605067_OBJECTS_H__