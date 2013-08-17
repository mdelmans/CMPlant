//
//  main.cpp
//  CMPlant
//
//  Created by Mihails Delmans on 25/07/2013.
//  Copyright (c) 2013 Mihails Delmans. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <OpenCL/opencl.h>
#include "CMPlant.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include "CLGeo.h"


#define N_IDX 0
#define DISP_H 800
#define DISP_W 1000

#define P_AREA  s[0]
#define P_THIC  s[1]
#define P_VISC  s[2]
#define P_YMOD  s[3]

#define P_SLEN  s[4]
#define P_SANG  s[5]
#define P_AEXK  s[6]
#define P_AVIS  s[7]



float area = 1, thick = 1, ymod = 1000, visc = 1000;

float pressure = 0.008, kEnv = 50, dt = 0.1;
float surfTen = 1;

cl_float8 newWallProperties = {area,thick,visc,ymod,0,surfTen,0,0};

CMPlant plant(pressure,kEnv,dt, newWallProperties);
bool play;
bool displayNodes = false;
bool displayForces = true;
bool displayNumbers = false;

int count = 1;
clock_t t = clock();


void delay(float sec){
    clock_t end_time = clock() + CLOCKS_PER_SEC * sec;
    while (clock() < end_time);
}


void display()
{
    plant.GLDisplay(DISP_W, DISP_H, displayNodes, displayNumbers, displayForces, 0.01);
    
}

/*! glut reshape callback function.  GLUT calls this function whenever
 the window is resized, including the first time it is created.
 You can use variables to keep track the current window size.
 */
void reshape(int width, int height)
{
    /* tell OpenGL we want to display in a recangle that is the
     same size as the window
     */
    glViewport(0,0,width,height);
    
    /* switch to the projection matrix */
    glMatrixMode(GL_PROJECTION);
    
    /* clear the projection matrix */
    glLoadIdentity();
    
    /* set the camera view, orthographic projection in 2D */
    gluOrtho2D(0,width,0,height);
    
    /* switch back to the model view matrix */
    glMatrixMode(GL_MODELVIEW);
}

void keyboard (unsigned char key, int x, int y){
    switch (key) {
        case 32:{
            int cCellsN = plant.cellsN;
            for (int i = 1; i <= cCellsN; i++)
                if (plant.cellDivide[i])
                    plant.Divide2(i);

            //plant.CLWriteStructure();
            glutPostRedisplay();
        }
            break;
        
        case 27:
            plant.CLStop();
            exit (0);
            break;
        case 13:
            play = !play;
            break;
        
        case 49:
            displayNodes = !displayNodes;
            glutPostRedisplay();
            break;
        case 50:
            displayNumbers = !displayNumbers;
             glutPostRedisplay();
            break;
        case 51:
            displayForces = !displayForces;
            glutPostRedisplay();
            break;
        case 112:
            plant.PrintAll();
            break;
        case 113:
            plant.pressure += 0.001;
            break;
        case 119:
            plant.pressure -= 0.001;
            break;
        case 97:
            std::cout << plant.pressure << " " << plant.CellArea(1) << "\n";
            break;
        default:
            break;
        
    }
}

void mousePassive(int x, int y){
    plant.GUICursorInCell(x-DISP_W/2, -y+DISP_H/2, 100);
    
    //std::cout << x-DISP_W/2 << " " <<  -y+DISP_H/2 << "\n";
    //std::cout << plant.GUICursorInCell(x-DISP_W/2, -y+DISP_H/2, 100) << "\n";
    
}

void mouseClick(int button, int state, int x, int y){
    
    static int cellIdx;
    
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN){
        cellIdx = plant.guiCellSelected;
    }
    
    if (button == GLUT_LEFT_BUTTON && state == GLUT_UP && cellIdx == plant.guiCellSelected && cellIdx!=-1){
        plant.Divide2(cellIdx);
    }
    
    if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN){
        cellIdx = plant.guiCellSelected;
    }
    
    if (button == GLUT_RIGHT_BUTTON && state == GLUT_UP && cellIdx == plant.guiCellSelected && cellIdx!=-1){
        plant.SetCellParams(cellIdx, 100000, 5000);
    }

    
}

void timer (int value){
    if (play) {
       
        plant.Step();
    }
    glutPostRedisplay();
    glutTimerFunc(10, timer, 0);
    
    for (int i = 1; i <= 6; i++) {
        std::cout << plant.GetAngle(i, 1) << " ";
    }
    std::cout << "\n";
}

int main(int argc, char ** argv)
{
    plant.CreateOneCell(6, 5000);
    
    //plant.CLCreateContext();
    
    //plant.ReadCellsFromFile("/Users/mdelmans/Desktop/CMPlant/XML/marchantia.xml", DISP_W, DISP_H, 100);
    /*
    std::cout << plant.cellsN << " " << plant.nodesN;
    plant.CLCreateContext();
    
    clock_t time = clock();
    
    for (int i = 1; i <= 100; i++){
        plant.CLCalculateForces(128);
        //plant.CLReadForce();
    }
    std::cout << (clock()-time) / (float)CLOCKS_PER_SEC << "\n";
    */
    /*
    plant.nodesN = 7;
    plant.cellsN = 1;
    
    plant.cellNodes[1][0] = 7;
    
    plant.cellNodes[1][1] = 1;
    plant.cellNodes[1][2] = 2;
    plant.cellNodes[1][3] = 3;
    plant.cellNodes[1][4] = 4;
    plant.cellNodes[1][5] = 5;
    plant.cellNodes[1][6] = 6;
    plant.cellNodes[1][7] = 7;

    
    plant.nodeCells[1][0] = {1,1};
    plant.nodeCells[1][1] = {1,1};
    
    
    plant.nodeCells[2][0] = {1,1};
    plant.nodeCells[2][1] = {1,2};
    
    
    plant.nodeCells[3][0] = {1,1};
    plant.nodeCells[3][1] = {1,3};
    
        
    plant.nodeCells[4][0] = {1,1};
    plant.nodeCells[4][1] = {1,4};
    
        
    plant.nodeCells[5][0] = {1,1};
    plant.nodeCells[5][1] = {1,5};
    
    
    plant.nodeCells[6][0] = {1,1};
    plant.nodeCells[6][1] = {1,6};
    
    plant.nodeCells[7][0] = {1,1};
    plant.nodeCells[7][1] = {1,7};
    
    

    
    plant.nodesPos[1] = {4141.75,-24663.1,0};
    plant.nodesPos[2] = {-3740.69,-27228.2,0};
    plant.nodesPos[3] = {-10770.8,-29515.8,0};
    plant.nodesPos[4] = {-9204.42,-38469.9,0};
    plant.nodesPos[5] = {-447.295,-36537.2,0};
    plant.nodesPos[6] = {6578.76,-34867.1,0};
    plant.nodesPos[7] = {5551.35,-30565.3,0};
   
    
    float area = 1, thick = 1, ymod = 1, visc = 1;
    
    plant.wallProperties[1][1] = {area,thick,visc,ymod, SecondNorm(plant.nodesPos[1] - plant.nodesPos[2]) ,0,0,0};
    plant.wallProperties[2][1] = {area,thick,visc,ymod, SecondNorm(plant.nodesPos[2] - plant.nodesPos[3]) ,0,0,0};
    plant.wallProperties[3][1] = {area,thick,visc,ymod, SecondNorm(plant.nodesPos[3] - plant.nodesPos[4]) ,0,0,0};
    plant.wallProperties[4][1] = {area,thick,visc,ymod, SecondNorm(plant.nodesPos[4] - plant.nodesPos[5]) ,0,0,0};
    plant.wallProperties[5][1] = {area,thick,visc,ymod, SecondNorm(plant.nodesPos[5] - plant.nodesPos[6]) ,0,0,0};
    plant.wallProperties[6][1] = {area,thick,visc,ymod, SecondNorm(plant.nodesPos[6] - plant.nodesPos[1]) ,0,0,0};
    
    plant.PrintAll();
    plant.Divide2(1);
    plant.PrintAll();
    */
    
    // Open GL
    // ######################################################################################################################################################
    
       
    glutInitWindowSize(DISP_W, DISP_H);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInit(&argc, argv);
    
    glutCreateWindow("CMPlant");
    
    glutDisplayFunc(display);
    
    glutReshapeFunc(reshape);
    
    glutKeyboardFunc(keyboard);
    
    glutPassiveMotionFunc(mousePassive);
    
    glutMouseFunc(mouseClick);
    
    glutTimerFunc(10, timer, 0);
    
    glClearColor(1,1,1,1);
    
    glutMainLoop();
    
    /*
    plant.nodesN = 7;
    plant.cellsN = 2;
    
    plant.cellNodes[1][0] = 4;
    plant.cellNodes[1][1] = 1;
    plant.cellNodes[1][2] = 2;
    plant.cellNodes[1][3] = 3;
    plant.cellNodes[1][4] = 4;
    
    plant.cellNodes[2][0] = 5;
    plant.cellNodes[2][1] = 4;
    plant.cellNodes[2][2] = 3;
    plant.cellNodes[2][3] = 5;
    plant.cellNodes[2][4] = 6;
    plant.cellNodes[2][5] = 7;
    
    plant.nodeCells[1][0] = {1,1};
    plant.nodeCells[1][1] = {1,1};
    
    plant.nodeCells[2][0] = {1,1};
    plant.nodeCells[2][1] = {1,2};
    
    plant.nodeCells[3][0] = {2,2};
    plant.nodeCells[3][1] = {1,3};
    plant.nodeCells[3][2] = {2,2};
    
    plant.nodeCells[4][0] = {2,2};
    plant.nodeCells[4][1] = {1,4};
    plant.nodeCells[4][2] = {2,1};
    
    plant.nodeCells[5][0] = {1,1};
    plant.nodeCells[5][1] = {2,3};
    
    plant.nodeCells[6][0] = {1,1};
    plant.nodeCells[6][1] = {2,4};
    
    plant.nodeCells[7][0] = {1,1};
    plant.nodeCells[7][1] = {2,5};
    
    plant.PrintAll();
    plant.InsertNewWall(1, 3, 1, 0, 0);
    plant.PrintAll();
    */
    
    
    return 0;
}



