///-lglut -lGLU -lGL
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <filesystem>
#include <GL/glut.h>

#include "geometry.h"
#include "objects.h"
#include "bitmap_image.hpp"


namespace GL
{
    ///Constants
    const double Delta = 10, AngleDelta = Geometry::PI / 30;
    const Geometry::Point POS(170, 0, 100), LOOK(-1, 0, 0);
    const Geometry::Point UP(0, 0, 1), RIGHT(0, 1, 0);
    const double fovy = 90, aspect_ratio = 1.33, near = 1, far = 10000;
    const double floorWidth = 1000, tileWidth = 20;
    const double windowWidth = 800, windowHeight = 600;
    const std::string outputPath="out/";

    const int MAX_OBJ = 1000;
    const Objects::Color BLACK{0, 0, 0}, WHITE{1, 1, 1};

    ///Data
    Geometry::Point pos, up, right, look;
    bool showaxis;
    int recursion_depth;

    int imageWidth, imageHeight;

    ///Keyboard/Mouse interaction
    void keyboardListener(unsigned char key, int x, int y);
    void specialKeyListener(int key, int x, int y);
    void mouseListener(int button, int state, int x, int y);
    
    ///Drawing
    void drawObjects();
    void drawAxis();
    void display();
    void animate();

    ///Initialize
    void loadData(std::string filename);
    void init(std::string filename);

    ///Capture
    void capture();

    ///Logging
    void printInfo();
    void printCameraInfo();
    void printAmbientInfo();
    void printSpecularInfo();
    void printDiffuseInfo();


    void printCameraInfo()
    {
        std::cout << "Camera Info: ";
        std::cout << "pos: " << pos << " ";
        std::cout << "look: " << look << " ";
        std::cout << "up: " << up << " ";
        std::cout << "right: " << right << " ";
        std::cout << std::endl;
    }

    void printAmbientInfo() {
        std::cout<<"Ambient Lighting "<<(Objects::showAmbient ? "On" : "Off")<<std::endl;
    }


    void printDiffuseInfo() {
        std::cout<<"Diffuse Lighting "<<(Objects::showDiffuse ? "On" : "Off")<<std::endl;
    }


    void printSpecularInfo() {
        std::cout<<"Specular Lighting "<<(Objects::showSpecular ? "On" : "Off")<<std::endl;
    }
    void printSpecularInfo();
    void printDiffuseInfo();
    void printInfo() {
        printCameraInfo();
        printAmbientInfo();
        printSpecularInfo();
        printDiffuseInfo();
    }

    //https://stackoverflow.com/a/10467633
    const std::string currentDateTime() {
        time_t     now = time(0);
        struct tm  tstruct;
        char       buf[80];
        tstruct = *localtime(&now);
        strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
        return buf;
    }

    Objects::Color getGLPixelColor(int x, int y){
        float pixel[3];
        glReadPixels(x, windowHeight-1-y, 1, 1, GL_RGB, GL_FLOAT, &pixel);
        return Objects::Color{pixel[0], pixel[1], pixel[2]};
    }

    void keyboardListener(unsigned char key, int x, int y)
    {            
        std::string s;
        using Geometry::PI;
        switch (key)
        {
        case '1':
            look = rotate(look, up, AngleDelta);
            right = rotate(right, up, AngleDelta);
            break;
        case '2':
            look = rotate(look, up, -AngleDelta);
            right = rotate(right, up, -AngleDelta);
            break;
        case '3':
            look = rotate(look, right, AngleDelta);
            up = rotate(up, right, AngleDelta);
            break;
        case '4':
            look = rotate(look, right, -AngleDelta);
            up = rotate(up, right, -AngleDelta);
            break;
        case '5':
            up = rotate(up, look, AngleDelta);
            right = rotate(right, look, AngleDelta);
            break;
        case '6':
            up = rotate(up, look, -AngleDelta);
            right = rotate(right, look, -AngleDelta);
            break;
        case 'x':
            pos = POS;
            look = LOOK;
            up = UP;
            right = RIGHT;
            break;
        case 'a':
            Objects::showAmbient ^= 1;
            printAmbientInfo();
            break;
        case 's':
            Objects::showSpecular ^= 1;
            printSpecularInfo();
            break;
        case 'd':
            Objects::showDiffuse ^= 1;
            printDiffuseInfo();
            break;
        case 'z':
            exit(0);
            break;
        case ' ':
            std::cout<<"Input Camera Info: ";
            std::cin>>s>>pos>>s>>look>>s>>up>>s>>right;
            break;
        case 'c':
        case '0':
            capture();
            break;
        default:
            break;
        }
    }

    void specialKeyListener(int key, int x, int y)
    {
        switch (key)
        {
        case GLUT_KEY_DOWN: //down arrow key
            pos = pos - look * Delta;
            break;
        case GLUT_KEY_UP: // up arrow key
            pos = pos + look * Delta;
            break;

        case GLUT_KEY_LEFT: //down arrow key
            pos = pos - right * Delta;
            break;
        case GLUT_KEY_RIGHT: // up arrow key
            pos = pos + right * Delta;
            break;

        case GLUT_KEY_PAGE_UP: //down arrow key
            pos = pos + up * Delta;
            break;
        case GLUT_KEY_PAGE_DOWN: // up arrow key
            pos = pos - up * Delta;
            break;
        default:
            break;
        }
    }

    void mouseListener(int button, int state, int x, int y)
    { //x, y is the x-y of the screen (2D)
        switch (button)
        {
        case GLUT_LEFT_BUTTON:
            if (state == GLUT_DOWN) {
                capture();
            }
            break;
        case GLUT_RIGHT_BUTTON:
            if (state == GLUT_DOWN)
                showaxis = !showaxis;
            break;
        case GLUT_MIDDLE_BUTTON:
            if (state == GLUT_DOWN)
                printInfo();
            break;
        default:
            break;
        }
    }

    void loadData(std::string filename)
    {
        std::ifstream sceneFile(filename, std::ifstream::in);
        assert(sceneFile);
        sceneFile >> recursion_depth >> imageWidth;
        imageHeight = imageWidth;

        int nObjects;
        sceneFile >> nObjects;
        assert(0 <= nObjects && nObjects <= MAX_OBJ);
        Objects::objects.resize(nObjects);

        for (int i = 0; i < nObjects; i++)
        {

            std::string type;
            sceneFile >> type;

            Objects::Object *obj = nullptr;
            std::cout << "Found " << type << std::endl;

            if (type == "sphere")
            {
                Geometry::Point center;
                double radius;
                sceneFile >> center >> radius;
                obj = new Objects::Sphere(center, radius);
            }
            else if (type == "triangle")
            {
                Geometry::Point A, B, C;
                sceneFile >> A >> B >> C;
                obj = new Objects::Triangle(A, B, C);
            }
            else if (type == "general")
            {
                double A, B, C, D, E, F, G, H, I, J;
                Geometry::Point ref;
                double length, height, width;
                sceneFile >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;
                sceneFile >> ref >> length >> width >> height;

                obj = new Objects::GeneralQuadricSurface(A, B, C, D, E, F, G, H, I, J, ref, length, width, height);
            }

            sceneFile >> obj->color;

            Objects::Reflection_coefficients cof;
            sceneFile >> cof.ambient >> cof.diffuse >> cof.specular >> cof.recursive >> cof.exponent;
            obj->coefficients = cof;

            if (obj != nullptr) Objects::objects[i] = obj;
        }
        int nLights;
        sceneFile >> nLights;
        assert(0 <= nLights && nLights <= MAX_OBJ);
        Objects::lights.resize(nLights);

        for (int i = 0; i < nLights; i++)
        {
            Objects::Light light;
            sceneFile >> light.position >> light.color;
            Objects::lights[i] = light;
        }
        Objects::objects.push_back(new Objects::Floor(floorWidth, tileWidth, BLACK, WHITE));
        sceneFile.close();
    }

    void init(std::string filename)
    {
        loadData(filename);
        pos = POS;
        look = LOOK;
        up = UP;
        right = RIGHT;
        showaxis = false;
        std::filesystem::create_directories("out/Past/");


        glClearColor(0, 0, 0, 0);     //clear the screen
        glMatrixMode(GL_PROJECTION); //load the PROJECTION matrix
        glLoadIdentity();             //initialize the matrix

        //fovy, aspect_ratio, near, far
        gluPerspective(fovy, aspect_ratio, near, far);
    }

	void clear() {
		for (auto x: Objects::objects)	delete x;
		Objects::objects.clear();
		Objects::lights.clear();
	}

    void drawObjects()
    {
        for (auto o : Objects::objects)    o->draw();
        for (auto l : Objects::lights)    l.draw();
    }

    void drawAxis()
    {
        glColor3f(1, 0, 0);
        glLineWidth(5);
        glBegin(GL_LINES);
        {
            glVertex3f(0, 0, 1000);
            glVertex3f(0, 0, -1000);
            glVertex3f(0, 1000, 1);
            glVertex3f(0, -1000, 1);
            glVertex3f(1000, 0, 1);
            glVertex3f(-1000, 0, 1);
        }
        glEnd();
    }

    void display()
    {
        //clear the display
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(0, 0, 0, 0); //color black
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glMatrixMode(GL_MODELVIEW); //load the correct matrix -- MODEL-VIEW matrix
        glLoadIdentity();            //initialize the matrix

        ///Camera Position, Direction, Up Direction
        Geometry::Point lookAt = pos + look;
        gluLookAt(pos.x, pos.y, pos.z, lookAt.x, lookAt.y, lookAt.z, up.x, up.y, up.z);
        glMatrixMode(GL_MODELVIEW); //again select MODEL-VIEW

        drawObjects();
        if (showaxis)
            drawAxis();

        //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
        glutSwapBuffers();
    }
    void animate()
    {
        //codes for any changes in Models, Camera
        glutPostRedisplay();
    }

    void capture()
    {
        using Geometry::Point, Objects::INF, Objects::Ray, Objects::Color;
        std::cout << "Capture: ";
        printCameraInfo();

        bitmap_image image(imageWidth, imageHeight);
        image.clear();
        
        double planeDistance  = (windowHeight/2.0) / tan(Geometry::degreeToRadian(fovy/2));
        double du = windowWidth/imageWidth;
        double dv = windowHeight/imageHeight;

        Point topleft = pos + look*planeDistance - right*windowWidth/2 + up*windowHeight/2;
        topleft = topleft + right*du/2 - up*dv/2;

        for (int x=0; x<imageWidth; x++) {
            for (int y=0; y<imageHeight; y++) {
                Point curPixel = topleft + right*du*x - up*dv*y;

                

                Ray ray(pos, curPixel-pos);
    

                double minDistance = INF;
                Color closestColor = BLACK;

                for (auto obj: Objects::objects) {
                    auto [distance, color] = obj -> intersect(ray, recursion_depth);
                    if (distance >= 0 && distance < minDistance) {
                        minDistance = distance;
                        closestColor = color;
                    }
                }

                //closestColor = getGLPixelColor(x*du, y*dv);
                image.set_pixel(x, y, closestColor.r*255, closestColor.g*255, closestColor.b*255);
            }
        }
        std::cout<<"Captured. Writing to File"<<std::endl;
        image.save_image(outputPath+"/capture.bmp");
        image.save_image(outputPath+"/Past/"+currentDateTime()+".bmp");
        std::cout<<"Finished."<<std::endl;
    }
}

int main(int argc, char **argv)
{
    std::string filename;
    if (argc > 1)
    {
        filename = argv[1];
    }
    else
    {
        filename = "scene.txt";
    }

    using namespace GL;
    glutInit(&argc, argv);
    glutInitWindowSize(windowWidth, windowHeight);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); //DEPTH, Double buffer, RGB color

    glutCreateWindow("My OpenGL Program");

    init(filename);

    glEnable(GL_DEPTH_TEST);  //enable DEPTH Testing
    glutDisplayFunc(display); //display callback function
    glutIdleFunc(animate);      //what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop(); //The main loop of OpenGL
	clear();

    return 0;
}
