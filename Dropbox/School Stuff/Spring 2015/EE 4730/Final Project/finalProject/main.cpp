#include <unistd.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <array>
#include <fstream>
#include <string>
#include <sstream>
#include "GLUT/GLUT.h"
#include "Mesh.h"
#include "Iterators.h"

using namespace std;

struct treeNode {
    struct treeNode * parent, * firstChild, *secondChild;
    Point center;
    double radius;
};

class BDTree {
public:
    Mesh * mesh;
    treeNode * head;
    int levels;
    Point pos;
    BDTree(Mesh * m, int l);
    ~BDTree();
    void treeInit(treeNode * node, vector<Point> points, int splitAxis);
    void deleteNode(treeNode * node);
    bool display(int i = 0);
    void levelCheck(treeNode * node, int currentLevel, int destLevel);
    void displayNode(treeNode * node);
};

BDTree::BDTree(Mesh * m, int l) {
    mesh = m;
    levels = l;
    pos.v[0] = 0;
    pos.v[1] = 0;
    pos.v[2] = 0;
    head = new treeNode;
    //head->parent = head;
    vector<Point> points;
    for (int i = 0; i < mesh->numVertices(); i++) {
        points.push_back(mesh->indVertex(i)->point());
    }
    treeInit(head, points, 0);
}

BDTree::~BDTree() {
    deleteNode(head);
}

void BDTree::deleteNode(treeNode *node) {
    if (node->firstChild)
        deleteNode(node->firstChild);
    if (node->secondChild)
        deleteNode(node->secondChild);
    delete node;
}

void BDTree::treeInit(treeNode * node, vector<Point> points, int level) {
    if (points.size() < 2) {
        node->radius = 0;
        if (level <= levels) {
            node->firstChild = new treeNode;
            node->secondChild = new treeNode;
            node->firstChild->parent = node;
            node->secondChild->parent = node;
            treeInit(node->firstChild, points, level + 1);
            treeInit(node->secondChild, points, level + 1);
        }
    }
    else {
        Point x, y, z, v;
        vector<Point> points1, points2;
        double distance = 0, temp;
        x = points[0];
        y = x;
        z = x;
        for (int i = 0; i < points.size(); i++) {
            v = points[i];
            temp = (x.v[0] - v.v[0]) * (x.v[0] - v.v[0]);
            temp += (x.v[1] - v.v[1]) * (x.v[1] - v.v[1]);
            temp += (x.v[2] - v.v[2]) * (x.v[2] - v.v[2]);
            temp = sqrt(temp);
            if (temp > distance) {
                distance = temp;
                y = v;
            }
        }
        distance = 0;
        for (int i = 0; i < points.size(); i++) {
            v = points[i];
            temp = (y.v[0] - v.v[0]) * (y.v[0] - v.v[0]);
            temp += (y.v[1] - v.v[1]) * (y.v[1] - v.v[1]);
            temp += (y.v[2] - v.v[2]) * (y.v[2] - v.v[2]);
            temp = sqrt(temp);
            if (temp > distance) {
                distance = temp;
                z = v;
            }
        }
        
        node->center.v[0] = (z.v[0] + y.v[0]) /2;
        node->center.v[1] = (z.v[1] + y.v[1]) /2;
        node->center.v[2] = (z.v[2] + y.v[2]) /2;
        node->radius = distance / 2;
        
        int i = 0;
        while (i < points.size()){
            v = points[i];
            temp = (node->center.v[0] - v.v[0]) * (node->center.v[0] - v.v[0]);
            temp += (node->center.v[1] - v.v[1]) * (node->center.v[1] - v.v[1]);
            temp += (node->center.v[2] - v.v[2]) * (node->center.v[2] - v.v[2]);
            temp = sqrt(temp);
            if (temp > node->radius && abs(temp - node->radius) > 0.001) {
                node->center.v[0] += ((temp - node->radius) / distance) * (v.v[0] - node->center.v[0]);
                node->center.v[1] += ((temp - node->radius) / distance) * (v.v[1] - node->center.v[1]);
                node->center.v[2] += ((temp - node->radius) / distance) * (v.v[2] - node->center.v[2]);
                node->radius = (node->radius + temp) / 2;
                i = 0;
            }
            else
                i++;
        }
        for (int i = 0; i < points.size(); i++) {
            if (points[i].v[level % 3] < node->center.v[level % 3])
                points1.push_back(points[i]);
            else
                points2.push_back(points[i]);
        }
        
        if (level <= levels) {
            node->firstChild = new treeNode;
            node->secondChild = new treeNode;
            node->firstChild->parent = node;
            node->secondChild->parent = node;
            treeInit(node->firstChild, points1, level + 1);
            treeInit(node->secondChild, points2, level + 1);
        }
    }
}

bool BDTree::display(int i) {
    if (i > levels) {
        cout << "Desired level is greated than number of levels.\n";
        return false;
    }
    else {
        treeNode * node = head;
        int level = 0;
        levelCheck(node, level, i);
        return true;
    }
}

void BDTree::levelCheck(treeNode * node, int curentLevel, int destLevel) {
    if (curentLevel == destLevel)
        displayNode(node);
    else {
        levelCheck(node->firstChild, curentLevel + 1, destLevel);
        levelCheck(node->secondChild, curentLevel + 1, destLevel);
    }
}

void BDTree::displayNode(treeNode * node) {
    if (node->radius == 0) {
        displayNode(node->parent);
    }
    else {
        glColor3f(0.5f, 0.5f, 0.7f);
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glTranslatef(node->center[0] + pos[0], node->center[1] + pos[1], node->center[2] + pos[2]);
        glutWireSphere(node->radius, 50, 50);
        glPopMatrix();
    }
}

BDTree * tree, * tree2;

array<float, 3> cameraPosition = {0,0,2}, objCenter = {0,0,0};
array<float, 3> boxMin = {0,0,0}, boxMax = {0,0,0};
vector<Point> * fNormal, * vNormal;
vector<double> * heTrgAngle;
float axislen = 1.414;

Mesh * pmesh, * pmesh2;
GLUquadricObj * obj;
float whratio;
int win_height, win_width;
int mousePositionX0 = 0, mousePositionY0 = 0;
int mouseButton = 0;
float my_Rotation_x = 0, my_Rotation_y = 0;
array<float, 3> my_Translation = { 0, 0, 0 };
array<double, 2> secPos = {.5, 0};
bool moveObj = 0;

void ComputeBoundingBox();
void ComputeNormal();

void glInit();
void display();
void reshape(int w, int h);
void mouseClick(int button, int state, int x, int y);
void mouseMove(int x, int y);
void keyboard(unsigned char key, int x, int y);
void setCamera();
void renderBoxesAndAxes();
void renderMesh();
void collision(treeNode * node1, treeNode * node2, int level);

int main(int argc, char** argv) {
    
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " input_mesh.mesh \n";
        return 1;
    }
    
    pmesh = new Mesh;
    pmesh2 = new Mesh;
    
    bool flag = pmesh->readMFile(argv[1]);
    if (!flag) {
        std::cerr << "Fail to read " << argv[1] << "\n";
        return 1;
    }
    
    pmesh->copyTo(* pmesh2);
    
    fNormal = new vector<Point>(pmesh->numFaces());
    vNormal = new vector<Point>(pmesh->numVertices());
    heTrgAngle = new vector<double>(pmesh->numEdges()*2);
    
    ComputeBoundingBox();
    ComputeNormal();
    
    for (int i = 0; i < pmesh2->numVertices(); i++)
        pmesh2->indVertex(i)->point().v[0] += 1;
    
    tree = new BDTree(pmesh, 10);
    tree2 = new BDTree(pmesh2, 10);
    
    
    glutInit(&argc, argv);
    glInit();
    
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouseClick);
    glutMotionFunc(mouseMove);
    glutKeyboardFunc(keyboard);
    
    glutMainLoop();
    
    delete fNormal;
    delete vNormal;
    delete heTrgAngle;
}
    
void ComputeBoundingBox() {
    objCenter[0]=objCenter[1]=objCenter[2]=0;
    boxMin[0]=boxMin[1]=boxMin[2]=1e5;
    boxMax[0]=boxMax[1]=boxMax[2]=-1e5;
    for (MeshVertexIterator vit(pmesh); !vit.end(); ++vit){
        Vertex * v = * vit;
        for (int j = 0; j < 3; ++j){
            float value = v->point()[j];
            objCenter[j]+= value;
            if (boxMax[j] < value)
                boxMax[j] = value;
            if (boxMin[j] > value)
                boxMin[j] = value;
        }
    }
    axislen=sqrt((boxMax[2]-boxMin[2])*(boxMax[2]-boxMin[2])+(boxMax[1]-boxMin[1])*(boxMax[1]-boxMin[1])+(boxMax[0]-boxMin[0])*(boxMax[0]-boxMin[0]));
    
    objCenter[0] /= pmesh->numVertices();
    objCenter[1] /= pmesh->numVertices();
    objCenter[2] /= pmesh->numVertices();
    
    cameraPosition[0]=objCenter[0];
    cameraPosition[1]=objCenter[1];
    cameraPosition[2]=objCenter[2]+axislen*1.5;
}

void ComputeNormal() {
    //First, compute face normals
    for (MeshFaceIterator fit(pmesh); !fit.end(); ++fit){
        Face * f = * fit;
        Halfedge * he = f->he();
        Halfedge * nhe = he->next();
        Vertex * v1 = he->source();
        Vertex * v2 = he->target();
        Vertex * v3 = nhe->target();
        Point nface1 = (v2->point() - v1->point());
        Point nface2 = (v3->point() - v2->point());
        Point nface = nface1 ^ nface2;
        nface /= nface.norm();
        (* fNormal)[f->index()] = nface;
    }
    
    //Second, compute corner angles
    //Need to set the index for all the halfedges; note that the index for faces and vertices is constructed when reading the mesh file
    for (MeshFaceIterator fit(pmesh); !fit.end(); ++fit){
        Face * f = *fit;
        array<Halfedge *, 3> he;
        array<double, 3> len, len2;
        he[0] = f->he();
        he[1] = he[0]->next();
        he[2] = he[0]->prev();
        for (int i = 0; i<3; ++i)
        {
            len2[i] = (he[i]->target()->point() - he[i]->source()->point()).norm2();
            len[i] = sqrt(len2[i]);
        }
        for (int i = 0; i<3; ++i)
        {
            int ip1 = (i + 1) % 3;
            int ip2 = (i + 2) % 3;
            if (len[i]<1e-8 || len[ip1]<1e-8)
                (* heTrgAngle)[he[i]->index()] = 3.14159265359 / 2;   //to handle the degenerate case
            else
                (* heTrgAngle)[he[i]->index()] = acos((len2[i] + len2[ip1] - len2[ip2]) / (2 * len[i] * len[ip1]));
        }
    }
    
    //Third, compute vertex normal
    for (MeshVertexIterator vit(pmesh); !vit.end(); ++vit)
    {
        Point sumVNorm(0, 0, 0);
        double sumWeight = 0;
        Vertex * v = * vit;
        for (VertexInHalfedgeIterator heit(v); !heit.end(); ++heit)
        {
            Halfedge * he = *heit;
            Face * f = he->face();
            int he_ind = he->index();
            int f_ind = f->index();
            
            Point temp = (* fNormal)[f_ind] * (* heTrgAngle)[he_ind];
            sumVNorm += temp;
            sumWeight += (* heTrgAngle)[he_ind];
        }
        int v_ind = v->index();
        (* vNormal)[v_ind] = sumVNorm / sumWeight;
    }
}

void glInit() {
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(900, 900);
    glutInitWindowPosition(50, 50);
    glutCreateWindow("Bounded Deformation Tree");
    glEnable(GL_CULL_FACE);
    glFrontFace(GL_CCW);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0,0,0,0);
    glShadeModel(GL_SMOOTH);
    glPolygonMode(GL_FRONT , GL_FILL);
    
    obj = gluNewQuadric();	//only for drawing spheres and cones
    
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    
    // Create light components
    GLfloat ambientLight0[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat	diffuseLight0[] = { 0.8f, 0.8f, 0.8, 1.0f };
    GLfloat specularLight0[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    GLfloat position0[] = {cameraPosition[0], cameraPosition[1], cameraPosition[2], 1.0f }; // the light is on the camera position
    
    // Assign created components to GL_LIGHT0
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight0);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight0);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight0);
    glLightfv(GL_LIGHT0, GL_POSITION, position0);
    
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
}

void display() {
    setCamera();
    //renderBoxesAndAxes();
    renderMesh();
    glutSwapBuffers();
}

void reshape(int w, int h) {
    glViewport (0, 0, (GLsizei) w, (GLsizei) h);
    win_height = h;
    win_width = w;
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
    whratio = (double)w / (double)h; 	//A Commonly Suggested Setting: set ratio in gluPerspective to the aspect ratio of the associated viewport
    gluPerspective(60, whratio, axislen*0.01, axislen*5);
    glMatrixMode (GL_MODELVIEW);	//change back to modelview
    glutPostRedisplay();

}

void mouseClick(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
        mouseButton = GLUT_LEFT_BUTTON;
    else if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN)
        mouseButton = GLUT_MIDDLE_BUTTON;
    else if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
        mouseButton = GLUT_RIGHT_BUTTON;
    mousePositionX0 = x;
    mousePositionY0 = y;
}

void mouseMove(int x, int y) {
    double movingScale = axislen / win_height;
    if (!moveObj) {
        /* rotation*/
        if (mouseButton == GLUT_LEFT_BUTTON ) {
            my_Rotation_y += x - mousePositionX0;
            my_Rotation_x += y - mousePositionY0;
        }
    
        /*xy translation */
        if (mouseButton == GLUT_MIDDLE_BUTTON) {
            my_Translation[0] += movingScale * (x - mousePositionX0);
            my_Translation[1] -= movingScale * (y - mousePositionY0);
        }
    
        /* zoom in and out */
        if (mouseButton == GLUT_RIGHT_BUTTON) {
            // suppose we want to make moving up as zooming out
            my_Translation[2] += movingScale * (y - mousePositionY0);
        }
    }
    else {
        for (int i = 0; i < pmesh2->numVertices(); i++) {
            pmesh2->indVertex(i)->point().v[0] += movingScale * (x - mousePositionX0);
            pmesh2->indVertex(i)->point().v[1] -= movingScale * (y - mousePositionY0);
        }
        tree2->pos.v[0] += movingScale * (x - mousePositionX0);
        tree2->pos.v[1] -= movingScale * (y - mousePositionY0);
    }
    mousePositionX0 = x;
    mousePositionY0 = y;
    glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
    if (key == 'm')
        moveObj = !moveObj;
}

void setCamera() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(cameraPosition[0], cameraPosition[1], cameraPosition[2], objCenter[0], objCenter[1], objCenter[2], 0, 1, 0);
    
    glTranslatef(my_Translation[0], my_Translation[1], my_Translation[2]);
    
    glTranslatef(objCenter[0], objCenter[1], objCenter[2]);	//before doing rotation to the object, move the object center to the origin
    
    glRotatef(my_Rotation_y, 0.0, 1.0, 0.0);
    glRotatef(my_Rotation_x, 1.0, 0.0, 0.0);
    
    glTranslatef(-objCenter[0], -objCenter[1], -objCenter[2]);
}

void renderBoxesAndAxes() {
    float axiswidth = axislen / 100;
    
    glMatrixMode (GL_MODELVIEW);
    
    glColor3f (1,1,1);
    
    glPushMatrix();
    //bounding box
    glTranslatef(objCenter[0], objCenter[1], objCenter[2]);
    glutWireCube(axislen);
    glTranslatef(-axislen/2, -axislen/2, -axislen/2);
    glutSolidSphere(axiswidth*1.5, 10, 10);
    
    //x-axis
    glColor3f (1,0,0);
    glPushMatrix();
    glRotatef(90, 0, 1, 0);
    gluCylinder(obj, axiswidth, axiswidth, axislen, 5, 5);
    glPopMatrix();
    
    glPushMatrix();
    glTranslatef(axislen,0,0);
    glRotatef(90, 0, 1, 0);
    glutWireCone(axiswidth*1.5, axiswidth*3, 10, 10);
    glPopMatrix();
    
    
    //y-axis
    glColor3f (0,1,0);
    glPushMatrix();
    glTranslatef(0,axislen,0);
    glRotatef(-90, 1, 0, 0);
    glutWireCone(axiswidth*1.5,axiswidth*3,10,10);
    glPopMatrix();
    
    glPushMatrix();
    glRotatef(-90, 1, 0, 0);
    gluCylinder(obj, axiswidth, axiswidth, axislen, 5, 5);
    glPopMatrix();
    
    //z-axis
    glColor3f (0,0,1);
    glPushMatrix();
    glTranslatef(0,0,axislen);
    glRotatef(-90, 0, 0, 1);
    glutWireCone(axiswidth*1.5,axiswidth*3,10,10);
    glPopMatrix();
    
    glPushMatrix();
    gluCylinder(obj, axiswidth, axiswidth, axislen, 5, 5);
    glPopMatrix();
    
    glPopMatrix();
}

void renderMesh() {
    glColor3f(0.7f, 0.7f, 0.7f);
    //traverse all the face and draw them
    glBegin(GL_TRIANGLES);
    for (MeshFaceIterator fit(pmesh); !fit.end(); ++fit){
        Face * f = *fit;
        Halfedge * he = f->he();
        array<Vertex *, 3> v;
        v[0] = he->source();
        v[1] = he->target();
        v[2] = he->next()->target();
        for (int i = 0; i < 3; ++i) {
            glNormal3dv((* vNormal)[v[i]->index()].v);
            glVertex3dv(v[i]->point().v);
        }
    }
    glEnd();
    glBegin(GL_TRIANGLES);
    for (MeshFaceIterator fit(pmesh2); !fit.end(); ++fit){
        Face * f = *fit;
        Halfedge * he = f->he();
        array<Vertex *, 3> v;
        v[0] = he->source();
        v[1] = he->target();
        v[2] = he->next()->target();
        for (int i = 0; i < 3; ++i) {
            glNormal3dv((* vNormal)[v[i]->index()].v);
            glVertex3dv(v[i]->point().v);
        }
    }
    glEnd();
    collision(tree->head, tree2->head, 0);
}

void collision(treeNode * node1, treeNode * node2, int level) {
    if (level < 5) {
        double d;
        d = (node1->center.v[0] - (node2->center.v[0] + tree2->pos[0]))
        * (node1->center.v[0] - (node2->center.v[0] + tree2->pos[0]));
        d += (node1->center.v[1] - (node2->center.v[1] + tree2->pos[1]))
        * (node1->center.v[1] - (node2->center.v[1] + tree2->pos[1]));
        d += (node1->center.v[2] - node2->center.v[2]) * (node1->center.v[2] - node2->center.v[2]);
        d = sqrt(d);
        if (d < node1->radius + node2->radius) {
            if (node1->firstChild && node2->firstChild)
                collision(node1->firstChild, node2->firstChild, level + 1);
            if (node1->firstChild && node2->secondChild)
                collision(node1->firstChild, node2->secondChild, level + 1);
            if (node1->secondChild && node2->firstChild)
                collision(node1->secondChild, node2->firstChild, level + 1);
            if (node1->secondChild && node2->secondChild)
                collision(node1->secondChild, node2->secondChild, level + 1);
        }
        else {
            if (node1 != tree->head || node2 != tree2->head) {
                tree->displayNode(node1->parent);
                tree2->displayNode(node2->parent);
            }
        }
            
    }
}