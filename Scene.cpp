#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include <iostream>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
    const char *str;
    XMLDocument xmlDoc;
    XMLElement *xmlElement;

    xmlDoc.LoadFile(xmlPath);

    XMLNode *rootNode = xmlDoc.FirstChild();

    // read background color
    xmlElement = rootNode->FirstChildElement("BackgroundColor");
    str = xmlElement->GetText();
    sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

    // read culling
    xmlElement = rootNode->FirstChildElement("Culling");
    if (xmlElement != NULL)
    {
        str = xmlElement->GetText();

        if (strcmp(str, "enabled") == 0)
        {
            this->cullingEnabled = true;
        }
        else
        {
            this->cullingEnabled = false;
        }
    }

    // read cameras
    xmlElement = rootNode->FirstChildElement("Cameras");
    XMLElement *camElement = xmlElement->FirstChildElement("Camera");
    XMLElement *camFieldElement;
    while (camElement != NULL)
    {
        Camera *camera = new Camera();

        camElement->QueryIntAttribute("id", &camera->cameraId);

        // read projection type
        str = camElement->Attribute("type");

        if (strcmp(str, "orthographic") == 0)
        {
            camera->projectionType = ORTOGRAPHIC_PROJECTION;
        }
        else
        {
            camera->projectionType = PERSPECTIVE_PROJECTION;
        }

        camFieldElement = camElement->FirstChildElement("Position");
        str = camFieldElement->GetText();
        sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

        camFieldElement = camElement->FirstChildElement("Gaze");
        str = camFieldElement->GetText();
        sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

        camFieldElement = camElement->FirstChildElement("Up");
        str = camFieldElement->GetText();
        sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

        camera->gaze = normalizeVec3(camera->gaze);
        camera->u = crossProductVec3(camera->gaze, camera->v);
        camera->u = normalizeVec3(camera->u);

        camera->w = inverseVec3(camera->gaze);
        camera->v = crossProductVec3(camera->u, camera->gaze);
        camera->v = normalizeVec3(camera->v);

        camFieldElement = camElement->FirstChildElement("ImagePlane");
        str = camFieldElement->GetText();
        sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
               &camera->left, &camera->right, &camera->bottom, &camera->top,
               &camera->near, &camera->far, &camera->horRes, &camera->verRes);

        camFieldElement = camElement->FirstChildElement("OutputName");
        str = camFieldElement->GetText();
        camera->outputFilename = string(str);

        this->cameras.push_back(camera);

        camElement = camElement->NextSiblingElement("Camera");
    }

    // read vertices
    xmlElement = rootNode->FirstChildElement("Vertices");
    XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
    int vertexId = 1;

    while (vertexElement != NULL)
    {
        Vec3 *vertex = new Vec3();
        Color *color = new Color();

        vertex->colorId = vertexId;

        str = vertexElement->Attribute("position");
        sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

        str = vertexElement->Attribute("color");
        sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

        this->vertices.push_back(vertex);
        this->colorsOfVertices.push_back(color);

        vertexElement = vertexElement->NextSiblingElement("Vertex");

        vertexId++;
    }

    // read translations
    xmlElement = rootNode->FirstChildElement("Translations");
    XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
    while (translationElement != NULL)
    {
        Translation *translation = new Translation();

        translationElement->QueryIntAttribute("id", &translation->translationId);

        str = translationElement->Attribute("value");
        sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

        this->translations.push_back(translation);

        translationElement = translationElement->NextSiblingElement("Translation");
    }

    // read scalings
    xmlElement = rootNode->FirstChildElement("Scalings");
    XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
    while (scalingElement != NULL)
    {
        Scaling *scaling = new Scaling();

        scalingElement->QueryIntAttribute("id", &scaling->scalingId);
        str = scalingElement->Attribute("value");
        sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

        this->scalings.push_back(scaling);

        scalingElement = scalingElement->NextSiblingElement("Scaling");
    }

    // read rotations
    xmlElement = rootNode->FirstChildElement("Rotations");
    XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
    while (rotationElement != NULL)
    {
        Rotation *rotation = new Rotation();

        rotationElement->QueryIntAttribute("id", &rotation->rotationId);
        str = rotationElement->Attribute("value");
        sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

        this->rotations.push_back(rotation);

        rotationElement = rotationElement->NextSiblingElement("Rotation");
    }

    // read meshes
    xmlElement = rootNode->FirstChildElement("Meshes");

    XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
    while (meshElement != NULL)
    {
        Mesh *mesh = new Mesh();

        meshElement->QueryIntAttribute("id", &mesh->meshId);

        // read projection type
        str = meshElement->Attribute("type");

        if (strcmp(str, "wireframe") == 0)
        {
            mesh->type = WIREFRAME_MESH;
        }
        else
        {
            mesh->type = SOLID_MESH;
        }

        // read mesh transformations
        XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
        XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

        while (meshTransformationElement != NULL)
        {
            char transformationType;
            int transformationId;

            str = meshTransformationElement->GetText();
            sscanf(str, "%c %d", &transformationType, &transformationId);

            mesh->transformationTypes.push_back(transformationType);
            mesh->transformationIds.push_back(transformationId);

            meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
        }

        mesh->numberOfTransformations = mesh->transformationIds.size();

        // read mesh faces
        char *row;
        char *cloneStr;
        int v1, v2, v3;
        XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
        str = meshFacesElement->GetText();
        cloneStr = strdup(str);

        row = strtok(cloneStr, "\n");
        while (row != NULL)
        {
            int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

            if (result != EOF)
            {
                mesh->triangles.push_back(Triangle(v1, v2, v3));
            }
            row = strtok(NULL, "\n");
        }
        mesh->numberOfTriangles = mesh->triangles.size();
        this->meshes.push_back(mesh);

        meshElement = meshElement->NextSiblingElement("Mesh");
    }
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
    this->image[i][j].r = c.r;
    this->image[i][j].g = c.g;
    this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
    if (this->image.empty())
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            vector<Color> rowOfColors;
            vector<double> rowOfDepths;

            for (int j = 0; j < camera->verRes; j++)
            {
                rowOfColors.push_back(this->backgroundColor);
                rowOfDepths.push_back(1.01);
            }

            this->image.push_back(rowOfColors);
            this->depth.push_back(rowOfDepths);
        }
    }
    else
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            for (int j = 0; j < camera->verRes; j++)
            {
                assignColorToPixel(i, j, this->backgroundColor);
                this->depth[i][j] = 1.01;
                this->depth[i][j] = 1.01;
                this->depth[i][j] = 1.01;
            }
        }
    }
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
    if (value >= 255.0)
        return 255;
    if (value <= 0.0)
        return 0;
    return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
    ofstream fout;

    fout.open(camera->outputFilename.c_str());

    fout << "P3" << endl;
    fout << "# " << camera->outputFilename << endl;
    fout << camera->horRes << " " << camera->verRes << endl;
    fout << "255" << endl;

    for (int j = camera->verRes - 1; j >= 0; j--)
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
                 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
                 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
        }
        fout << endl;
    }
    fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
    string command;

    // TODO: Change implementation if necessary.
    command = "./magick convert " + ppmFileName + " " + ppmFileName + ".png";
    system(command.c_str());
}

double Scene::absOfDouble(double a)
{
    if(a < 0)
    {
        return -a;
    }
    else
    {
        return a;
    }
}

// Finds smallest component of given three
double Scene::findSmallestComponent(double x, double y, double z)
{
    double absX = absOfDouble(x);
    double absY = absOfDouble(y);
    double absZ = absOfDouble(z);

    double min = absX;
    double resultMin = x;
    if(absY < min)
    {
        min = absY;
        resultMin = y;
    }
    if(absZ < min)
    {
        min = absZ;
        resultMin = z;
    }

    return resultMin;
}

Matrix4 Scene::findModelingTransformationMatrix(Mesh& currentMesh)
{
    // call this function for each mesh

    int numOfTransformations = currentMesh.numberOfTransformations;

    Matrix4 transformationMatrix = getIdentityMatrix();

    for(int j = 0; j < numOfTransformations; j++)
    {
        int currentTransformationId = currentMesh.transformationIds[j] - 1;
        char currentTransformationType = currentMesh.transformationTypes[j];

        if(currentTransformationType == 't')
        {
            double x = this->translations[currentTransformationId]->tx;
            double y = this->translations[currentTransformationId]->ty;
            double z = this->translations[currentTransformationId]->tz;

            double translationMatrix[4][4] = {
                    {1,0,0,x},
                    {0,1,0,y},
                    {0,0,1,z},
                    {0,0,0,1}
            };

            Matrix4 tMatrix(translationMatrix);


            transformationMatrix = multiplyMatrixWithMatrix(tMatrix,transformationMatrix);

        }
        else if(currentTransformationType == 's')
        {
            double x = this->scalings[currentTransformationId]->sx;
            double y = this->scalings[currentTransformationId]->sy;
            double z = this->scalings[currentTransformationId]->sz;

            double scalingMatrix[4][4] = {
                    {x,0,0,0},
                    {0,y,0,0},
                    {0,0,z,0},
                    {0,0,0,1}
            };

            Matrix4 sMatrix(scalingMatrix);

            transformationMatrix = multiplyMatrixWithMatrix(sMatrix,transformationMatrix);

        }
        else if(currentTransformationType == 'r')
        {
            double ang = this->rotations[currentTransformationId]->angle;
            double x = this->rotations[currentTransformationId]->ux;
            double y = this->rotations[currentTransformationId]->uy;
            double z = this->rotations[currentTransformationId]->uz;

            Vec3 u,v,w;

            u = Vec3(x,y,z);

            double smallestComp = findSmallestComponent(x,y,z);
            if(x == smallestComp)
            {
                v = Vec3(0,z,-y);
            }
            else if(y == smallestComp)
            {
                v = Vec3(z,0,-x);
            }
            else // z == smallestComp
            {
                v = Vec3(y,-x,0);
            }

            w = crossProductVec3(u,v);

            v = normalizeVec3(v);
            w = normalizeVec3(w);

            double tempRotationMatrix[4][4] =  {
                    {u.x,u.y,u.z,0},
                    {v.x,v.y,v.z,0},
                    {w.x,w.y,w.z,0},
                    {0,0,0,1}
            };

            Matrix4 tempRMatrix(tempRotationMatrix);

            double inverseTempRotationMatrix[4][4] =  {
                    {u.x,v.x,w.x,0},
                    {u.y,v.y,w.y,0},
                    {u.z,v.z,w.z,0},
                    {0,0,0,1}
            };

            Matrix4 inverseTempRMatrix(inverseTempRotationMatrix);

            double xRotationMatrix[4][4] =  {
                    {1,0,0,0},
                    {0,cos(ang * M_PI / 180),-sin(ang * M_PI / 180),0},
                    {0,sin(ang * M_PI / 180),cos(ang * M_PI / 180),0},
                    {0,0,0,1}
            };

            Matrix4 xRMatrix(xRotationMatrix);

            Matrix4 finalRotationMatrix = multiplyMatrixWithMatrix(xRMatrix,tempRMatrix);
            finalRotationMatrix = multiplyMatrixWithMatrix(inverseTempRMatrix,finalRotationMatrix);
            transformationMatrix = multiplyMatrixWithMatrix(finalRotationMatrix,transformationMatrix);
        }

    }


    return transformationMatrix;

}
Matrix4 Scene::CameraTransformationMatrix(Camera &camera)
{
    Matrix4 result;
    Vec3 cu = camera.u;
    Vec3 cv = camera.v;
    Vec3 cw = camera.w;
    Vec3 e = camera.position;
    result.values[0][0] = cu.x;
    result.values[0][1] = cu.y;
    result.values[0][2] = cu.z;
    result.values[0][3] = -(cu.x*e.x + cu.y*e.y + cu.z*e.z);

    result.values[1][0] = cv.x;
    result.values[1][1] = cv.y;
    result.values[1][2] = cv.z;
    result.values[1][3] = -(cv.x*e.x + cv.y*e.y + cv.z*e.z);

    result.values[2][0] = cw.x;
    result.values[2][1] = cw.y;
    result.values[2][2] = cw.z;
    result.values[2][3] = -(cw.x*e.x + cw.y*e.y + cw.z*e.z);

    result.values[3][3] = 1;

    return result;
}

Matrix4 Scene::P2O(Camera &camera)
{
    double f = camera.far, n = camera.near;
    Matrix4 result;
    result.values[0][0] = n;
    result.values[1][1] = n;
    result.values[2][2] = f+n;
    result.values[2][3] = f*n;
    result.values[3][2] = -1;

    return result;
}

Matrix4 Scene::OrthographicProjection(Camera &camera)
{
    double r = camera.right, l = camera.left, t = camera.top, b = camera.bottom, f = camera.far, n = camera.near;
    Matrix4 result;
    result.values[0][0] = 2/(r-l);
    result.values[0][3] = (r+l)/(l-r);

    result.values[1][1] = 2/(t-b);
    result.values[1][3] = (t+b)/(b-t);

    result.values[2][2] = 2/(n-f);
    result.values[2][3] = (f+n)/(n-f);

    result.values[3][3] = 1;

    return result;
}

void Scene::PerspectiveDivide(Vec4 &v)
{
    double w = v.t;
    v.x = v.x/w;
    v.y = v.y/w;
    v.z = v.z/w;
    v.t  = 1;
}

Matrix4 Scene::ViewportTransformation(int nx, int ny) // nx = camera->horRes, ny = camera->verRes
{
    Matrix4 result;
    result.values[0][0] = nx/2;
    result.values[0][3] = (nx-1)/2;

    result.values[1][1] = ny/2;
    result.values[1][3] = (ny-1)/2;

    result.values[2][2] = 0.5;
    result.values[2][3] = 0.5;

    result.values[3][3] = 1;

    return result;
}

Matrix4 Scene::findViewingTransformationsMatrix(Camera &camera)
{
    Matrix4 viewingMatrix = getIdentityMatrix();

//	double camMatrixMove[4][4] = {
//								{1,0,0,-camera.position.x},
//								{0,1,0,-camera.position.y},
//								{0,0,1,-camera.position.z},
//								{0,0,0,1}
//								};
//
//	double camMatrixRotate[4][4] = {
//								{camera.u.x,camera.u.y,camera.u.z,0},
//								{camera.v.x,camera.v.y,camera.v.z,0},
//								{camera.w.x,camera.w.y,camera.w.z,0},
//								{0,0,0,1}
//								};
//
//	Matrix4 cameraTransformationMatrix = multiplyMatrixWithMatrix(camMatrixRotate,camMatrixMove);
    Matrix4 cameraTransformationMatrix = CameraTransformationMatrix(camera);
    Matrix4 orthographicMatrix = OrthographicProjection(camera);

//	double ortographicMatrix[4][4] = {
//								{(2/(camera.right-camera.left)),0,0,-((camera.right+camera.left)/(camera.right-camera.left))},
//								{0,(2/(camera.top-camera.bottom)),0,-((camera.top+camera.bottom)/(camera.top-camera.bottom))},
//								{0,0,-(2/(camera.far-camera.near)),-((camera.far+camera.near)/(camera.far-camera.near))},
//								{0,0,0,1}
//								};

//	double p2oMatrix[4][4] = {
//							{camera.near,0,0,0},
//							{0,camera.near,0,0},
//							{0,0,(camera.far + camera.near),(camera.far * camera.near)},
//							{0,0,-1,0}
//							};





    if(camera.projectionType == 0) // ortographic
    {
        viewingMatrix = multiplyMatrixWithMatrix(cameraTransformationMatrix,viewingMatrix);
        viewingMatrix = multiplyMatrixWithMatrix(orthographicMatrix,viewingMatrix);
    }
    else if(camera.projectionType == 1) // perspective
    {
        Matrix4 p2oMatrix = P2O(camera);
        Matrix4 perspectiveMatrix = multiplyMatrixWithMatrix(orthographicMatrix,p2oMatrix);
        viewingMatrix = multiplyMatrixWithMatrix(cameraTransformationMatrix,viewingMatrix);
        viewingMatrix = multiplyMatrixWithMatrix(perspectiveMatrix,viewingMatrix);
    }

    return viewingMatrix;
}

Vec4 Scene::FindViewportCoordinates(Vec4 point, Camera &camera) // matrix here is the transformation matrix up to perspective divide.
{
    //PerspectiveDivide(point);
    return multiplyMatrixWithVec4(ViewportTransformation(camera.horRes, camera.verRes), point);
}

bool Scene::lbVisible(double den, double num, double &te, double &tl)
{
    double t;

    if(den > 0) // potentially entering
    {
        t = num / den;
        if(t > tl)
        {
            return false;
        }
        if(t > te)
        {
            te = t;
        }
    }
    else if(den < 0) // potentially leaving
    {
        t = num / den;
        if(t < te)
        {
            return false;
        }
        if(t < tl)
        {
            tl = t;
        }
    }
    else if(num > 0) // line parallel to edge
    {
        return false;
    }

    return true;
}

bool Scene::LiangBarskyAlgorithm(Vec4 &vertexZero, Vec4 &vertexOne, Color &c0, Color&c1, Camera &camera)
{
    double dx = vertexOne.x - vertexZero.x;
    double dy = vertexOne.y - vertexZero.y;
    double dz = vertexOne.z - vertexZero.z;

    double tE = 0;
    double tL = 1;
    bool visible = false;
//    Color c0 = *colorsOfVertices[vertexZero.colorId];
//    Color c1 = *colorsOfVertices[vertexOne.colorId];
    Color dc = c1 - c0;

    if(lbVisible(dx, (-1) - vertexZero.x, tE, tL) &&
       lbVisible(-dx, vertexZero.x - (1), tE, tL) &&
       lbVisible(dy, (-1) - vertexZero.y, tE, tL) &&
       lbVisible(-dy, vertexZero.y - (1), tE, tL) &&
       lbVisible(dz, (-1) - vertexZero.z, tE, tL) &&
       lbVisible(-dz, vertexZero.z - (1), tE, tL))
    {
        visible = true;
        if(tL < 1)
        {
            vertexOne.x = vertexZero.x + dx*tL;
            vertexOne.y = vertexZero.y + dy*tL;
            vertexOne.z = vertexZero.z + dz*tL;
            c1 = c0 + (tL * dc);
        }
        if(tE > 0)
        {
            vertexZero.x = vertexZero.x + dx*tE;
            vertexZero.y = vertexZero.y + dy*tE;
            vertexZero.z = vertexZero.z + dz*tE;
            c0 = c0 + (tE * dc);
        }
    }

    return visible;
}


bool Scene::isBackfaceCulling(Vec4& v1, Vec4& v2, Vec4& v3, Camera& camera)
{
    // check if culling is enabled !

    Vec3 v1_3(v1.x,v1.y,v1.z);
    Vec3 v2_3(v2.x,v2.y,v2.z);
    Vec3 v3_3(v3.x,v3.y,v3.z);

    // v2-v1
    Vec3 u = subtractVec3(v2_3,v1_3);

    //v3-v1
    Vec3 w = subtractVec3(v3_3,v1_3);

    Vec3 normalVec = normalizeVec3(crossProductVec3(u,w));

    Vec3 center( ((v1_3.x + v2_3.x + v3_3.x)/3) , ((v1_3.y + v2_3.y + v3_3.y)/3) , ((v1_3.z + v2_3.z + v3_3.z)/3) );

    //Vec3 v = normalizeVec3(subtractVec3(center,camera.position));

    double facingValue = dotProductVec3(normalVec,center);

    if(facingValue > 0)
    {
        return false; // back face culling is not applied (visible)
    }
    else
    {
        return true; // back face culling is applied (not visible)
    }
}

double Scene::ComputeLineEquations(Vec4& v1, Vec4& v2, double x, double y)
{
    double result = ( x*(v1.y - v2.y) ) + ( y*(v2.x - v1.x) ) + (v1.x*v2.y) - (v1.y*v2.x);

    //std::cout << result << endl;
    return result;
}


void Scene::TriangleRasterization(Vec4& v1, Vec4& v2, Vec4& v3, Camera& camera)
{
    

    double yMin = v1.y;
    double yMax = v1.y;
    double xMin = v1.x;
    double xMax = v1.x;

    // set yMin
    if(v2.y < yMin)
    {
        yMin = v2.y;
    }

    if(v3.y < yMin)
    {
        yMin = v3.y;
    }

    // set yMax
    if(v2.y > yMax)
    {
        yMax = v2.y;
    }

    if(v3.y > yMax)
    {
        yMax = v3.y;
    }

    // set xMin
    if(v2.x < xMin)
    {
        xMin = v2.x;
    }

    if(v3.x < xMin)
    {
        xMin = v3.x;
    }

    // set xMax
    if(v2.x > xMax)
    {
        xMax = v2.x;
    }

    if(v3.x > xMax)
    {
        xMax = v3.x;
    }

    // apply algorithm
    double alpha = 0;
    double beta = 0;
    double gamma = 0;
    Color currentColor;   

    // i: y values
    // j: x values
    for(int i = yMin; i <= yMax; i++)
    {
        for(int j = xMin; j <= xMax; j++)
        {
            alpha = ComputeLineEquations(v2,v3,j,i) / ComputeLineEquations(v2,v3,v1.x,v1.y);
            beta = ComputeLineEquations(v3,v1,j,i) / ComputeLineEquations(v3,v1,v2.x,v2.y);
            gamma = ComputeLineEquations(v1,v2,j,i) / ComputeLineEquations(v1,v2,v3.x,v3.y);

            if( alpha >= 0 && beta >= 0 && gamma >= 0)
            {
                currentColor.r = alpha * this->colorsOfVertices[v1.colorId]->r +
                                 beta * this->colorsOfVertices[v2.colorId]->r +
                                 gamma * this->colorsOfVertices[v3.colorId] ->r;

                currentColor.g = alpha * this->colorsOfVertices[v1.colorId]->g +
                                 beta * this->colorsOfVertices[v2.colorId]->g +
                                 gamma *this->colorsOfVertices[v3.colorId] ->g;

                currentColor.b = alpha * this->colorsOfVertices[v1.colorId]->b +
                                 beta * this->colorsOfVertices[v2.colorId]->b +
                                 gamma *this->colorsOfVertices[v3.colorId] ->b;

                double totalDepth = alpha*v1.z + beta*v2.z + gamma*v3.z;

                DepthBufferTest(j,i,currentColor,totalDepth,camera);
            }
        }
    }

}


void Scene::Midpoint(Vec4& v0,Vec4& v1, Color c0, Color c1, Camera &camera) {
    double dx = v1.x - v0.x;
    double dy = v1.y - v0.y;
    double z, z0 = v0.z, z1 = v1.z, dz;
    int d, slopeSign = 1;
    Color dc, c;

    if (abs(dy) <= abs(dx)) {
        if (v1.x < v0.x) {
            swap(v0, v1);
            swap(c0, c1);
            swap(z0,z1);
        }
        if (v1.y < v0.y)
            slopeSign = -1;

        int y = (int) v0.y;
        c = c0;
        z = z0;
        d = 2*(v0.y - v1.y) + (slopeSign * (v1.x - v0.x));
        dc = (c1 - c0) / (v1.x - v0.x);
        dz = (z1-z0) / (v1.x - v0.x);
        for (int x = v0.x; x < v1.x; x++) {
            DepthBufferTest(x,y,c,z,camera);
            d += 2*(v0.y - v1.y); // choose E or NE, this is gonna happen
            if (d * slopeSign < 0) { // choose NE
                y += slopeSign;
                d += 2*(slopeSign * (v1.x - v0.x));
            }

            c = c + dc;
            z += dz;
        }
    }
    else if (abs(dy) > abs(dx)) {
        if (v1.y < v0.y) {
            swap(v0, v1);
            swap(c0, c1);
            swap(z0,z1);
        }
        if (v1.x < v0.x) {
            slopeSign = -1;
        }

        int x = v0.x;
        c = c0;
        z = z0;
        d = 2*(v1.x - v0.x) + (slopeSign * (v0.y - v1.y));
        dc = (c1 - c0) / (v1.y - v0.y);
        dz = (z1-z0) / (v1.y - v0.y);
        for (int y = v0.y; y < v1.y; y++) {
            DepthBufferTest(x,y,c,z,camera);
            d += 2*(v1.x - v0.x);
            if (d * slopeSign > 0) {
                x += slopeSign;
                d += 2*(slopeSign * (v0.y - v1.y));
            }
            c = c + dc;
            z += dz;
        }
    }
}

// TODO: depth buffer algorithm
void Scene::DepthBufferTest(int x, int y, Color currentColor, double currentDepth, Camera& camera)
{
    if(x < 0 || x >= camera.horRes || y < 0 || y >= camera.verRes)
    {
        return;
    }

    double bufferDepth = this->depth[x][y];

    currentColor.r = makeBetweenZeroAnd255(round(currentColor.r));
    currentColor.g = makeBetweenZeroAnd255(round(currentColor.g));
    currentColor.b = makeBetweenZeroAnd255(round(currentColor.b));

    

    // overwrite the buffer
    if(currentDepth <= bufferDepth)
    {
        
        this->depth[x][y] = currentDepth;
        
        //std::cout << currentColor.r << endl;

        // paint the color matrix
        assignColorToPixel(x,y,currentColor);
        
    }
}

/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
    // TODO: Implement this function

    Matrix4 viewingMatrix = findViewingTransformationsMatrix(*camera);

    

    for(auto currentMesh: this->meshes)
    {
        // default backface culling: disabled
        bool isBackfaceCullingEnabled = false;
        bool isClippingApplicable = false;

        // set backface culling mode
        if(this->cullingEnabled)
        {
            isBackfaceCullingEnabled = true;
        }

        // wireframe mesh: do clipping
        if(currentMesh->type == 0)
        {
            isClippingApplicable = true;
        }

        Matrix4 meshModelingMatrix = findModelingTransformationMatrix(*currentMesh);


        Matrix4 meshModelViewProjectionMatrix = multiplyMatrixWithMatrix(viewingMatrix,meshModelingMatrix);

        
        for(auto currentTriangle: currentMesh->triangles)
        {
            Vec3* vertex1 = this->vertices[currentTriangle.vertexIds[0] - 1];
            Vec4 v1(vertex1->x,vertex1->y,vertex1->z, 1, currentTriangle.vertexIds[0] - 1);

            Vec3* vertex2 = this->vertices[currentTriangle.vertexIds[1] - 1];
            Vec4 v2(vertex2->x,vertex2->y,vertex2->z, 1, currentTriangle.vertexIds[1] - 1);

            Vec3* vertex3 = this->vertices[currentTriangle.vertexIds[2] - 1];
            Vec4 v3(vertex3->x,vertex3->y,vertex3->z, 1, currentTriangle.vertexIds[2] - 1);
            //std::cout << currentMesh->meshId << endl;
            

            v1 = multiplyMatrixWithVec4(meshModelViewProjectionMatrix,v1);
            v2 = multiplyMatrixWithVec4(meshModelViewProjectionMatrix,v2);
            v3 = multiplyMatrixWithVec4(meshModelViewProjectionMatrix,v3);


            // check back face culling
            if(isBackfaceCullingEnabled)
            {
                if(isBackfaceCulling(v1,v2,v3,*camera))
                {
                    continue;
                }
            }

            // wireframe mode
            if(isClippingApplicable)
            {
                //pers divide
                PerspectiveDivide(v1);
                PerspectiveDivide(v2);
                PerspectiveDivide(v3);

                Vec4 v1Copy1(v1.x,v1.y,v1.z,v1.t,v1.colorId);
                Vec4 v1Copy2(v1.x,v1.y,v1.z,v1.t,v1.colorId);
                Vec4 v2Copy1(v2.x,v2.y,v2.z,v2.t,v2.colorId);
                Vec4 v2Copy2(v2.x,v2.y,v2.z,v2.t,v2.colorId);
                Vec4 v3Copy1(v3.x,v3.y,v3.z,v3.t,v3.colorId);
                Vec4 v3Copy2(v3.x,v3.y,v3.z,v3.t,v3.colorId);
                Color c1 = *colorsOfVertices[v1.colorId];
                Color c2 = *colorsOfVertices[v2.colorId];
                Color c3 = *colorsOfVertices[v3.colorId];
                // Line1: v1-v2
                bool isLine1Visible = LiangBarskyAlgorithm(v1Copy1,v2Copy1, c1, c2, *camera);

                // Line2: v2-v3
                bool isLine2Visible = LiangBarskyAlgorithm(v2Copy2,v3Copy1, c2 ,c3, *camera);

                // Line3: v3-v1
                bool isLine3Visible = LiangBarskyAlgorithm(v3Copy2,v1Copy2, c3 ,c1, *camera);


                // without perspective divide and viewport transform
                //v1 = FindViewportCoordinates(v1,*camera);
                //v2 = FindViewportCoordinates(v2,*camera);
                //v3 = FindViewportCoordinates(v3,*camera);

                
                v1Copy1 = FindViewportCoordinates(v1Copy1,*camera);
                v2Copy1 = FindViewportCoordinates(v2Copy1,*camera);
                v3Copy1 = FindViewportCoordinates(v3Copy1,*camera);

                v1Copy2 = FindViewportCoordinates(v1Copy2,*camera);
                v2Copy2 = FindViewportCoordinates(v2Copy2,*camera);
                v3Copy2 = FindViewportCoordinates(v3Copy2,*camera);
                


                if(isLine1Visible)
                {
                    Midpoint(v1Copy1,v2Copy1, c1 , c2, *camera);
                }

                if(isLine2Visible)
                {
                    Midpoint(v2Copy2,v3Copy1,c2, c3,*camera);
                }

                if(isLine3Visible)
                {
                    Midpoint(v3Copy2,v1Copy2,c3, c1,*camera);
                }
            }
            else // solid mode
            {
                // perspective divide and viewport transform

                //pers divide
                PerspectiveDivide(v1);
                PerspectiveDivide(v2);
                PerspectiveDivide(v3);
                
                v1 = FindViewportCoordinates(v1,*camera);
                v2 = FindViewportCoordinates(v2,*camera);
                v3 = FindViewportCoordinates(v3,*camera);

            
                TriangleRasterization(v1,v2,v3,*camera);
            }
        }

    }

    // for each model find modeling transformation matrix
    // multiply modeling transformation matrix with camera transformation matrix
    // for each triangle and vertex multiply with final matrix


}
