//
// Created by Eloi on 12/19/20.
//

#ifndef DGTCM_SOLVER_H
#define DGTCM_SOLVER_H


#include "PhysicalModel.h"
#include "Matrix.h"


namespace DGTCM
{
    class Solver {
        enum enum_direction {
            north = 1,
            south = 2,
            east = 3,
            west = 4,
            top = 5,
            bottom = 6
        };

    public:
        Solver(double itrEpsilon, double itrTimeDelta, double endTempDiff);

        virtual ~Solver();
        void Initialize();
        void Solve();
        void setTimeDelta(double deltaT) { m_timedelta = deltaT; };
        void setEpsilon(double epsilon) { m_epsilon = epsilon; };
        void Step();
    private:
        // Solver functions

        void ComputeTn1(double & oMaxTempInPiece);
        void ComputeQn();

        // Solver auxiliary functions
        double ComputeAiTi(unsigned int x, unsigned int y, unsigned int z);
        double ComputeAp(unsigned int x, unsigned int y, unsigned int z);
        double ComputeBp(unsigned int x, unsigned int y, unsigned int z);

        // Boundary conditions
        double ApBoundaryCondition(unsigned int x, unsigned int y, unsigned int z);
        double BpBoundaryCondition(unsigned int x, unsigned int y, unsigned int z);
        double QnBoundaryCondition(unsigned int x, unsigned int y, unsigned int z);

        // Helper functions
        Matrix* ActiveMatrix() { return *m_pActiveTemp; };
        double A(unsigned int x, unsigned int y, unsigned int z, int direction);
        double A(unsigned int x, unsigned int y, unsigned int z, enum_direction direction);
        double Cp(unsigned int x,  unsigned int y,  unsigned int z);
        double Alpha(unsigned int x,  unsigned int y,  unsigned int z);
        double Lambda(unsigned int x,  unsigned int y,  unsigned int z);
        double tempAt(unsigned int x,  unsigned int y,  unsigned int z);
        double tempAtDirection(unsigned int x,  unsigned int y,  unsigned int z, int direction);
        double tempAtDirection(unsigned int x,  unsigned int y,  unsigned int z, enum_direction direction);
        void SetActiveTempMatrix(Matrix *& iMatrix) { m_pActiveTemp = &iMatrix; };
        double GetSurfaceArea(unsigned int x, unsigned int y, unsigned int z);

// Data members
        double m_endTempDiff;
        double m_epsilon;
        double m_timedelta;
        Matrix** m_pActiveTemp;
        Matrix* m_Tn; // Temperature at time
        Matrix* m_Tsup; // Temperature at time
        Matrix* m_Tn1; // Temperature at time + delta time
        Matrix* m_Qn; // Net heat flow at state n
        const PhysicalModel* m_physical;
    };
}




#endif //DGTCM_SOLVER_H
