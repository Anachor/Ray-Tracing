#ifndef __1605067_OBJECTS_H__
#define __1605067_OBJECTS_H__

#include "1605067_geometry.h"
#include <tuple>
#include <algorithm>

namespace Objects
{
    const double INF = 1e100;
    bool showAmbient  = true;
    bool showSpecular = true;
    bool showDiffuse  = true;
    using namespace Geometry;

    struct Color;
    struct Ray;
    struct Reflection_coefficients;
    struct Object;
    struct Light;

    std::vector<Object*> objects;
    std::vector<Light> lights;

    struct Color
    {
        double r, g, b; //Between [0, 1]
        Color(double r, double g, double b) : r(r), g(g), b(b) {}

        Color() {}
        Color operator+(const Color &u) const
        {
            return Color(r + u.r, g + u.g, b + u.b);
        }
        Color operator-(const Color &u) const
        {
            return Color(r - u.r, g - u.g, b - u.b);
        }
        Color operator*(const double u) const
        {
            return Color(r * u, g * u, b * u);
        }
        Color operator/(const double u) const
        {
            return Color(r / u, g / u, b / u);
        }
        Color operator*(const Color &u) const
        {
            return Color(r * u.r, g * u.r, b * u.r);
        }
        friend std::istream &operator>>(std::istream &is, Color &o)
        {
            return is >> o.r >> o.g >> o.b;
        }
        friend std::ostream &operator<<(std::ostream &os, const Color &o)
        {
            return os << o.r << " " << o.g << " " << o.b;
        }

        void clamp() {
            r = std::clamp(r, 0.0, 1.0);
            g = std::clamp(g, 0.0, 1.0);
            b = std::clamp(b, 0.0, 1.0);
        }
    };

    struct Ray
    {
        Point a, dir;
        Ray(Point a, Point dir) : a(a), dir(Geometry::unit(dir)) {}
    };

    struct Reflection_coefficients
    {
        double ambient, diffuse, specular, recursive, exponent;
    };

    struct Light
    {

        Point position;
        Color color;

        void draw()
        {
            glPointSize(6.0);
            glColor3f(color.r, color.g, color.b);
            glBegin(GL_POINTS);
            {
                glVertex3f(position.x, position.y, position.z);
            }
            glEnd();
        }
    };

    struct Object
    {
        Color color;
        Reflection_coefficients coefficients;

        Object(Color color, Reflection_coefficients coefficients)
            : color(color), coefficients(coefficients)
        {}

        Object() {}

        virtual bool onSurface(Point p) = 0;
        virtual Point normal(Point p) = 0;
        virtual std::pair<double, Color> intersect(Ray r) = 0;
        virtual void draw() = 0;
        virtual ~Object() {}

        std::pair<double, Color> intersect(Ray cameraRay, int level)
        {
            if (level == 0)
                return this->intersect(cameraRay);

            auto [distance, basecolor] = this->intersect(cameraRay);
            if (distance == INF)        return  {INF, {}};
            
            Point X = cameraRay.a + Geometry::unit(cameraRay.dir)*distance; 
            
            Color color(0, 0, 0);
            
            ///Ambient
            if (showAmbient) color =  color + basecolor * coefficients.ambient;

            for (auto l: lights) {
                Ray lightRay(l.position, X-l.position);
                double lightDis = Geometry::length(X-l.position);
                bool blocked = false;

                Point N = this->normal(X);
                Point L = Geometry::unit(l.position - X);
                Point V = Geometry::unit(cameraRay.a - X);
                Point R = N*2*Geometry::dot(L,N) - L;

                
                if (Geometry::dot(L,N) < 0)     N=N*-1;
                if (Geometry::dcmp(Geometry::dot(V,N)) <= 0)     continue;

                for (auto obj: objects) {
                    auto [nwdis, color] = obj -> intersect(lightRay);
                    if (nwdis >= 0 && Geometry::dcmp(nwdis-lightDis) < 0) {
                        blocked = true;
                        break;
                    }
                }

                if (!blocked) {

                    ///Diffuse
                    double lambert = Geometry::dot(L, N);
                    if (showDiffuse) 
                        color = color + l.color*coefficients.diffuse*lambert*basecolor;

                    ///Specular

                    double phong = Geometry::dot(R, V);
                    if (phong > 0 && showSpecular) 
                        color = color + l.color * coefficients.specular * pow(phong, coefficients.exponent)*basecolor;
                }
            }

            color.clamp();
            return {distance, color};
        };
    };

    struct Floor : Object
    {
        double floorWidth;
        double tileWidth;
        Color darkcolor, lightcolor;

        Floor(double floorWidth, double tileWidth, Color darkcolor, Color lightcolor)
            : floorWidth(floorWidth), tileWidth(tileWidth), darkcolor(darkcolor), lightcolor(lightcolor)
        {
            coefficients = {0.2, 0.2, 0.2, 0.2, 20};
        }
        

        bool onSurface(Point P) override {
            return Geometry::dcmp(P.z) == 0 &&
                   abs(P.x) <= floorWidth / 2 &&
                   abs(P.y) <= floorWidth / 2;
        }
        Point normal(Point p) override {
            return Point(0, 0, 1);
        }

        std::pair<double, Color> intersect(Ray r) override
        {
            double d = -r.a.z / r.dir.z;

            Point P = r.a + r.dir * d;
            if (!onSurface(P))  return {INF, {}};
            
            int i = (P.x + floorWidth / 2) / tileWidth;
            int j = (P.y + floorWidth / 2) / tileWidth;

            std::pair<double, Color> ans;
            if (i % 2 == j % 2)
                ans = {d, darkcolor};
            else
                ans = {d, lightcolor};
            return ans;
        };

        void draw() override
        {
            assert(floorWidth > 0 && tileWidth > 0 && floorWidth / tileWidth <= 1000);
            for (auto [x, i] = std::make_tuple(-floorWidth / 2, 0); x <= floorWidth / 2; x += tileWidth, i++)
            {
                for (auto [y, j] = std::make_tuple(-floorWidth / 2, 0); y <= floorWidth / 2; y += tileWidth, j++)
                {
                    Color c = (i % 2 == j % 2 ? darkcolor : lightcolor);
                    glColor3f(c.r, c.g, c.b);

                    glBegin(GL_QUADS);
                    {
                        glVertex3f(x, y, 0);
                        glVertex3f(x + tileWidth, y, 0);
                        glVertex3f(x + tileWidth, y + tileWidth, 0);
                        glVertex3f(x, y + tileWidth, 0);
                    }
                    glEnd();
                }
            }
        }
    };

    struct Sphere : Object
    {
        Point center;
        double radius;

        const static int stacks = 100, slices = 100;

        Sphere(Point center, double radius) : center(center), radius(radius)
        {}

        bool onSurface(Point P) override {
            return dcmp(Geometry::length(P-center)) == 0;
        }
        Point normal(Point p) override {
            return Geometry::unit(p-center);
        }


        std::pair<double, Color> intersect(Ray r) override
        {
            // std::cout<<"Ray: "<<r.a<<" --> "<<r.dir<<"  ----  Sphere: O = "<<center<<", r = "<<radius<<std::endl;

            Geometry::Line L = {r.a, r.dir};
            double OX = Geometry::distancePointLine(center, L);
            double OA = length(r.a - center);

            double D = sqrt(OA * OA - OX * OX);
            double d = sqrt(radius * radius - OX * OX);

            Point Z1 = r.a + unit(r.dir) * D;
            Point Z2 = r.a - unit(r.dir) * D;
            if (Geometry::length(Z2 - center) < Geometry::length(Z1 - center))
                D = -D;

            std::pair<double, Color> ans;
            if (OX > radius)
                return {INF, {}};
            else if (D + d < 0)
                return {INF, {}};
            else if (D - d > 0)
                ans = {D - d, color};
            else
                ans = {D + d, color};
            return ans;
        }

        void draw() override
        {
            // write codes for drawing sphere
            using Geometry::Point, Geometry::PI;
            Point points[stacks + 1][slices + 1];

            //generate points
            for (int i = 0; i <= stacks; i++)
            {
                double h = radius * sin(((double)i / (double)stacks) * (PI / 2));
                double r = radius * cos(((double)i / (double)stacks) * (PI / 2));
                for (int j = 0; j <= slices; j++)
                {
                    points[i][j].x = r * cos(((double)j / (double)slices) * 2 * PI);
                    points[i][j].y = r * sin(((double)j / (double)slices) * 2 * PI);
                    points[i][j].z = h;
                }
            }
            //draw quads using generated points
            for (int i = 0; i < stacks; i++)
            {
                for (int j = 0; j < slices; j++)
                {
                    glBegin(GL_QUADS);
                    {
                        glColor3f(color.r, color.g, color.b);

                        //upper hemisphere
                        glVertex3f(center.x + points[i][j].x, center.y + points[i][j].y, center.z + points[i][j].z);
                        glVertex3f(center.x + points[i][j + 1].x, center.y + points[i][j + 1].y, center.z + points[i][j + 1].z);
                        glVertex3f(center.x + points[i + 1][j + 1].x, center.y + points[i + 1][j + 1].y, center.z + points[i + 1][j + 1].z);
                        glVertex3f(center.x + points[i + 1][j].x, center.y + points[i + 1][j].y, center.z + points[i + 1][j].z);
                        //lower hemisphere
                        glVertex3f(center.x + points[i][j].x, center.y + points[i][j].y, center.z - points[i][j].z);
                        glVertex3f(center.x + points[i][j + 1].x, center.y + points[i][j + 1].y, center.z - points[i][j + 1].z);
                        glVertex3f(center.x + points[i + 1][j + 1].x, center.y + points[i + 1][j + 1].y, center.z - points[i + 1][j + 1].z);
                        glVertex3f(center.x + points[i + 1][j].x, center.y + points[i + 1][j].y, center.z - points[i + 1][j].z);
                    }
                    glEnd();
                }
            }
        }
    };

    struct Triangle : Object
    {
        Point A, B, C;
        Plane plane;

        Triangle(Point A, Point B, Point C) : A(A), B(B), C(C), plane(A, B, C)
        {}

        bool onSurface(Point O) override {
            if (!plane.onPlane(O))  return false;

            double ABC = length(cross(A - B, B - C));
            double OAB = length(cross(A - O, B - O));
            double OBC = length(cross(B - O, C - O));
            double OCA = length(cross(C - O, A - O));
            double rat = (OAB+OBC+OCA)/ABC;
            
            return Geometry::dcmp(rat-1) == 0;
        }

        Point normal(Point p) override {
            return unit(plane.normal);
        }

        std::pair<double, Color> intersect(Ray r) override
        {
            using namespace Geometry;
            std::pair<double, Color> ans;

            double d = plane.intersectLine({r.a, r.dir}) / length(r.dir);
            if (d <= 0)      return {INF, {}};

            Point O = r.a + unit(r.dir) * d;
            if (onSurface(O))
                ans = {d, color};
            else
                return {INF, {}};

            return ans;
        };
        void draw() override
        {
            // write codes for drawing triangle
            glColor3f(color.r, color.g, color.b);
            glBegin(GL_TRIANGLES);
            {
                glVertex3f(A.x, A.y, A.z);
                glVertex3f(B.x, B.y, B.z);
                glVertex3f(C.x, C.y, C.z);
            }
            glEnd();
        }
    };

    std::vector<double> solveQuadratic(double A, double B, double C)
    {
        using Geometry::EPS;
        if (abs(A) < EPS)
        {
            if (abs(B) < EPS)
                return {};
            return {-C / B};
        }
        else
        {
            double D = B * B - 4 * A * C;
            if (D < -EPS)
                return {};
            else if (abs(D) < EPS)
                return {-B / (2 * A)};
            else
                return {(-B - sqrt(D)) / (2 * A), (-B + sqrt(D)) / (2 * A)};
        }
    }

    struct GeneralQuadricSurface : Object
    {
        ///Ax^2 + By^2 + Cz^2 + Dxy + Exz + Fyz + Gx + Hy + Iz + J
        double A, B, C, D, E, F, G, H, I, J;
        Geometry::Point reference_point;
        double length, width, height;

        GeneralQuadricSurface(double A, double B, double C, double D, double E,
                              double F, double G, double H, double I, double J,
                              Geometry::Point reference_point,
                              double length, double width, double height)
            : A(A), B(B), C(C), D(D), E(E), F(F), G(G), H(H), I(I), J(J),
              reference_point(reference_point), length(length), width(width), height(height)
        {}
        

        bool inRange(double x, double len, double X)
        {
            return len == 0 || (x <= X && X <= x + len);
        }

        bool inCube(Point Z)
        {
            return inRange(reference_point.x, length, Z.x) &&
                   inRange(reference_point.y, width, Z.y) &&
                   inRange(reference_point.z, height, Z.z);
        }

        bool onSurface(Point P)
        {
            double z = A * P.x * P.x + B * P.y * P.y + C * P.z * P.z + 
                       D * P.x * P.y + E * P.x * P.z + F * P.y * P.z + 
                       G * P.x + H * P.y + I * P.z + J;
            return Geometry::dcmp(z) == 0;
        }
        Point normal(Point P)
        {
            double ax = 2*A*P.x + D*P.y + E*P.z + G;
            double ay = 2*B*P.y + D*P.x + F*P.z + H;
            double az = 2*C*P.z + E*P.x + F*P.y + I;
            Point ans(ax, ay, az);
            return Geometry::unit(ans);
        }

        std::pair<double, Color> intersect(Ray r) override
        {

            double qA = A * r.dir.x * r.dir.x + B * r.dir.y * r.dir.y + C * r.dir.z * r.dir.z + 
                        D * r.dir.x * r.dir.y + E * r.dir.x * r.dir.z + F * r.dir.y * r.dir.z;

            double qB = 2 * A * r.a.x * r.dir.x + 2 * B * r.a.y * r.dir.y + 2 * C * r.a.z * r.dir.z + 
                        D * (r.a.x * r.dir.y + r.a.y * r.dir.x) + 
                        E * (r.a.z * r.dir.x + r.a.x * r.dir.z) + 
                        F * (r.a.y * r.dir.z + r.a.z * r.dir.y) + 
                        G * r.dir.x + H * r.dir.y + I * r.dir.z;

            double qC = A * r.a.x * r.a.x + B * r.a.y * r.a.y + C * r.a.z * r.a.z + 
                        D * r.a.x * r.a.y + E * r.a.x * r.a.z + F * r.a.y * r.a.z + 
                        G * r.a.x + H * r.a.y + I * r.a.z + J;

            double closest = INF, len = Geometry::length(r.dir);
            for (double d : solveQuadratic(qA, qB, qC)) {
                Point p = r.a + r.dir * d;
                if (!inCube(p))     continue;
                d = d * len;
                if (d < closest && d >= 0) {
                    closest = d;
                }
            }

            std::pair<double, Color> ans = {closest, color};
            if (closest == INF) return ans;

            return ans;
        };
        void draw() override
        {
            //std::cout<<"Cannot draw general quadratic surface"<<std::endl;
        }
    };

}

#endif // __1605067_OBJECTS_H__