//-----------------------------------------------------------------------
// Header file of SistemSPH class which simulate Vicsek model using
// Smoothed-particle hydrodynamics (SPH) method.
//
// Licensing: This code is distributed under the Apache License 2.0
// Author: Carlos Planelles Alemany, planelles20(at)gmail(dot)com
//-----------------------------------------------------------------------

#ifndef SYSTEM_SPH
#define SYSTEM_SPH

#include <iostream>
#include <string>
#include <fstream>
// GLEW
#define GLEW_STATIC
#include <GL/glew.h>

#include <cuda_gl_interop.h>
// cuda rand
#include<curand.h>
#include<curand_kernel.h>

//structs
struct particle {
    // particle index
    GLint id;
    //cell index
    GLint cellIdx;
    //position
    GLfloat x;
    GLfloat y;
    GLfloat z;
    //velocity
    GLfloat vx;
    GLfloat vy;
    GLfloat vz;
    //angles
    GLfloat theta;
    GLfloat alpha;

};

struct intvec2 {
    GLint init;
    GLint end;
};

struct calcStruct {
    GLfloat thetaMedium;
    GLfloat alphaMedium;
    GLint N; // # of neighbors
};

struct point {
    GLfloat x;
    GLfloat y;
    GLfloat z;
};

// sistem SPH class
class SystemSPH {
    public:
        SystemSPH(unsigned int blocks,
                  unsigned int threads,
                  unsigned int xMesh,
                  unsigned int yMesh,
                  unsigned int zMesh);

        virtual ~SystemSPH();
        void Particle_print();
        void Indices_print();
        void Calc_print();
        void Save(const std::string& nameFile);
        void Calculate();
        // draw methods (openGL)
        void DrawParticles();
        void BackGround(float r, float g, float b, float a);
        void DrawBoundary();
        void SeedUpdate(int i);

    protected:

    private:
        //atributes
        unsigned int numBlocks, numThreads;
        unsigned int numParticles;
        unsigned int xMeshDim, yMeshDim, zMeshDim;
        unsigned int numIndices;

        // buffers to simulation
        particle *h_particle; //host psticles
        intvec2 *h_inidices; //host indices
        calcStruct *h_calc; //host culculation

        particle *d_particle; //device psticles
        intvec2 *d_inidices; //device indices //NOTE: just to bugs
        calcStruct *d_calc; //device calculation

        // buffers to openGL
        //particles adn indiceCells buffers
        GLuint VAO, EBO;
        GLuint VBOparticles;
        GLuint VBOindices;
        GLuint VBOcalc;   //buffer for calculation
        cudaGraphicsResource_t cudaResourceBufParticles;
        cudaGraphicsResource_t cudaResourceBufIndices;
        cudaGraphicsResource_t cudaResourceBufCalc;
        GLushort *particleIndices;

        // Boundary
        GLuint VAOboundary, VBOboundary, EBOboundary;
        point *boundaryPoints;
        GLushort *boundaryIndices;

        ///////////////// Seedm to generate random numbers //////////////////
        GLfloat seed;

        //calculate methods
        void InitParticleData();
        void PosParticleCell();
        void SortParticles();
        void ClearGridIndices();
        void BuiltGridIncices();
        void CreateGridIndices();
        void CalcOperations();
        //plot methods
        void PolygonMode();
        void CreateBounderiesPoints();
        void CreateBounderiesIndices();
};

#endif // SYSTEM_SPH
