//
// Created by Eloi on 12/19/20.
//

#include "Solver.h"
#include <iostream>

using namespace DGTCM;

Solver::Solver(double itrEpsilon, double itrTimeDelta, double endTempDiff):
m_pActiveTemp(nullptr),
m_Tn(nullptr),
m_Tsup(nullptr),
m_Tn1(nullptr),
m_Qn(nullptr),
m_physical(nullptr),
m_epsilon(itrEpsilon),
m_timedelta(itrTimeDelta),
m_endTempDiff(endTempDiff){

}

Solver::~Solver() {
    if(m_Tn) {
        delete m_Tn;
        m_Tn = nullptr;
    }
    if(m_Tn1) {
        delete m_Tn1;
        m_Tn1 = nullptr;
    }
    if(m_Qn) {
        delete m_Qn;
        m_Qn = nullptr;
    }
    if(m_physical) {
        delete m_physical;
        m_physical = nullptr;
    }

}
void Solver::Initialize() {
    m_physical = new PhysicalModel(0.25, 0.2, 0.1, 50, 40, 20, 1200.0, 300.0, 7833.0, 1);
    m_Tn = new Matrix(m_physical->nX(), m_physical->nY(), m_physical->nZ(), m_physical->initialTemp());
    m_Tn1 = new Matrix(m_physical->nX(), m_physical->nY(), m_physical->nZ(), m_physical->initialTemp());
    m_Tsup = new Matrix(m_physical->nX(), m_physical->nY(), m_physical->nZ(), m_physical->initialTemp());
    m_Qn = new Matrix(m_physical->nX(), m_physical->nY(), m_physical->nZ(), 0.0);
}

void Solver::Solve() {
    Initialize();
    double maxTempInPiece = m_physical->initialTemp();
    double totalTime = 0.0;
    double maxTime = 200.0;
    while(totalTime < maxTime && maxTempInPiece > (m_physical->fluidTemp()+30)) {
        totalTime += m_timedelta;
        ComputeQn();
        ComputeTn1(maxTempInPiece);
        std::cout<<"Current Elapsed time: "<<totalTime<<std::endl;
        std::cout<<"Current Max internal temperature: "<<maxTempInPiece<<std::endl;
        m_Tn->Copy(m_Tn1);
    }
}


void Solver::ComputeTn1(double & oMaxTempInPiece) {
    m_Tsup->Copy(m_Tn1);
    SetActiveTempMatrix(m_Tsup);
    double maxDiff = 0.0;
    oMaxTempInPiece = 0.0;
    for(unsigned int x = 0; x < m_Tn1->XSize(); x++) {
        for (unsigned int y = 0; y < m_Tn1->YSize(); y++) {
            for (unsigned int z = 0; z < m_Tn1->ZSize(); z++) {
                // T * Ap = sum(lamda *T) + bp = AiTi * bp
                double AiTi = ComputeAiTi(x,y,z);
                double Ap = ComputeAp(x,y,z);
                double Bp = ComputeBp(x,y,z);
                m_Tn1->operator()(x,y,z) = (AiTi + Bp) / Ap;
                // Calculate max difference
                double val = abs(m_Tn1->at(x,y,z) - m_Tsup->at(x,y,z));
                if(val > maxDiff)
                    maxDiff = val;
                if(m_Tn1->at(x,y,z) > oMaxTempInPiece) oMaxTempInPiece = m_Tn1->at(x,y,z);
            }
        }
    }
    if (maxDiff > m_epsilon) return ComputeTn1(oMaxTempInPiece);
}


void Solver::ComputeQn() {
    SetActiveTempMatrix(m_Tn);
    // Do raw sum of a_side * T_side + boundary conditions
    for(unsigned int x = 0; x < m_Tn->XSize(); x++) {
        for (unsigned int y = 0; y < m_Tn->YSize(); y++) {
            for (unsigned int z = 0; z < m_Tn->ZSize(); z++) {
                double Q = 0.0;
                for (int direction = 1; direction <= 6; direction++) {
                    Q += A(x, y, z, direction)  * (tempAtDirection(x, y, z, direction) - tempAt(x,y,z));
                }
                Q += QnBoundaryCondition(x, y, z);
                m_Qn->operator()(x, y, z) = Q;
            }
        }
    }

}

double Solver::ComputeAiTi(unsigned int x, unsigned int y, unsigned int z) {
    double AiTi = 0.0;
    for (int direction = 1; direction <= 6; direction++) {
        AiTi += A(x, y, z, direction) * tempAtDirection(x, y, z,direction);
    }
    return m_physical->beta()*AiTi;
}

double Solver::ComputeAp(unsigned int x, unsigned int y, unsigned int z) {
    double Ap = m_physical->rho()*m_physical->vol()*Cp(x,y,z) / m_timedelta;
    for (int direction = 1; direction <= 6; direction++) {
        Ap += m_physical->beta() * A(x, y, z, direction);
    }
    Ap += ApBoundaryCondition(x,y,z);
    return Ap;
}

double Solver::ComputeBp(unsigned int x, unsigned int y, unsigned int z) {
    double Bp = (1-m_physical->beta())*m_Qn->operator()(x,y,z) + m_Tn->operator()(x,y,z) * m_physical->rho()*m_physical->vol()*Cp(x,y,z) / m_timedelta;

    Bp += BpBoundaryCondition(x,y,z);
    return Bp;
}

// Heat trtansfer at boundary = alpha * A Tn-Tfluid
double Solver::QnBoundaryCondition(unsigned int x, unsigned int y, unsigned int z) {
    return Alpha(x,y,z)*GetSurfaceArea(x,y,z)*(m_physical->fluidTemp() -  m_Tn->operator()(x,y,z));
}

double Solver::ApBoundaryCondition(unsigned int x, unsigned int y, unsigned int z) {
    return m_physical->beta()*Alpha(x,y,z)*GetSurfaceArea(x,y,z);
}

double Solver::BpBoundaryCondition(unsigned int x, unsigned int y, unsigned int z) {
    return m_physical->beta()*Alpha(x,y,z)*m_physical->fluidTemp()*GetSurfaceArea(x,y,z);
}


double Solver::tempAt(unsigned int x, unsigned int y, unsigned int z) {
    if (x <0 || x>ActiveMatrix()->XSize()-1) return 0.0;
    if (y <0 || y>ActiveMatrix()->YSize()-1) return 0.0;
    if (z <0 || z>ActiveMatrix()->ZSize()-1) return 0.0;
    return ActiveMatrix()->operator()(x,y,z);
}


double Solver::tempAtDirection(unsigned int x, unsigned int y, unsigned int z, int direction) {
    return tempAtDirection(x,y,z,enum_direction(direction));
}

double Solver::tempAtDirection(unsigned int x, unsigned int y, unsigned int z, enum_direction direction) {
    double T = 0.0;
    switch (direction) {
        case north:
            if (x < m_physical->nX()-1) T = tempAt(x+1, y, z);
            break;
        case south:
            if (x>0)                    T = tempAt(x-1, y, z);
            break;
        case east:
            if (y < m_physical->nY()-1) T = tempAt(x, y+1, z);
            break;
        case west:
            if (y>0)                    T = tempAt(x, y-1, z);
            break;
        case top:
            if (z < m_physical->nZ()-1) T = tempAt(x, y, z+1);
            break;
        case bottom:
            if (z>0)                    T = tempAt(x, y, z-1);
            break;
    }
    return T;
}

double Solver::A(unsigned int x, unsigned int y, unsigned int z, int direction) {
    return A(x, y, z, enum_direction(direction));
}

double Solver::A(unsigned int x, unsigned int y, unsigned int z, enum_direction direction) {
    double A = 0.0;
    switch (direction) {
        case north:
            if (x < m_physical->nX()-1) A = ((Lambda(x+1, y, z) + Lambda(x, y, z)) *  m_physical->xArea()) / (2.0 * m_physical->dX());
            break;
        case south:
            if (x>0)                    A = ((Lambda(x-1, y, z) + Lambda(x, y, z)) *  m_physical->xArea()) / (2.0 * m_physical->dX());
            break;
        case east:
            if (y < m_physical->nY()-1) A = ((Lambda(x , y+1, z) + Lambda(x, y, z)) *  m_physical->yArea()) / (2.0 * m_physical->dY());
            break;
        case west:
            if (y>0)                    A = ((Lambda(x , y-1, z) + Lambda(x, y, z)) *  m_physical->yArea()) / (2.0 * m_physical->dY());
            break;
        case top:
            if (z < m_physical->nZ()-1) A = ((Lambda(x , y, z+1) + Lambda(x, y, z)) *  m_physical->zArea()) / (2.0 * m_physical->dZ());
            break;
        case bottom:
            if (z>0)                    A = ((Lambda(x , y, z-1) + Lambda(x, y, z)) *  m_physical->zArea()) / (2.0 * m_physical->dZ());
            break;
    }
    return A;
}


double Solver::GetSurfaceArea(unsigned int x, unsigned int y, unsigned int z) {
    double Area = 0.0;
    if (x == 0 || x == (m_Tn->XSize()-1))
        Area += m_physical->xArea();
    if (y == 0 || y == (m_Tn->YSize()-1))
        Area += m_physical->yArea();
    if (z == 0 || z == (m_Tn->ZSize()-1))
        Area += m_physical->zArea();
    return Area;
}


double Solver::Alpha(unsigned int x, unsigned int y, unsigned int z) {
    //Alpha = -5561.455 + 53.22458*temp - 0.1940006*temp^2 + 0.0003336*temp^3 - 2.548565e-7*temp^4 + 6.988594e-11*temp^5
    double temp = ActiveMatrix()->operator()(x,y,z);
    /*
     if(temp>1273.0) { return 350.0;}
    if(temp>1173.15 && temp<=1273.15) { return 1100.0;}
    if(temp>1073.15 && temp<=1173.15) { return 1800.0;}
    if(temp>973.15 && temp<=1073.15) { return 2400.0;}
    if(temp>873.15  && temp<=973.15) { return 2600.0;}
    if(temp>673.15  && temp<=873.15)  { return 2300.0;}
    if(temp>573.15  && temp<=673.15)  { return 1500.0;}
    if(temp>473.15  && temp<=573.15)  { return 900.0;}
    if(temp>373.15  && temp<=473.15)  { return 550.0;}
    if(temp>273.15  && temp<=373.15)  { return 150.0;}
      */
    //return 400.0;
    double a = -5561.455 + 53.22458*temp - 0.1940006*pow(temp,2.0) + 0.0003336*pow(temp,3.0) - 2.548565e-7*pow(temp,4.0) + 6.988594e-11*pow(temp,5.0);
    return a;
}

double Solver::Cp(unsigned int x, unsigned int y, unsigned int z) {
    double temp = tempAt(x,y,z);
    double cp = 487.6;// - 1.57*temp + 2*pow(10.0, -2.0)*pow(temp, 2.0) - 7.69 * pow(10.0, -5.0) * pow(temp, 3.0) +
    //1.2*pow(10.0, -7.0)*pow(temp, 4.0) - 6.67 * pow(10.0,-11.0)*pow(temp, 5.0);
    return cp;
}


double Solver::Lambda(unsigned int x, unsigned int y, unsigned int z) {
    double l = 29.0;
    double temp = ActiveMatrix()->operator()(x,y,z);
    if (temp < 1273.0) {
        l += (1273.0-temp)/900.0 * 29;
    }
    return l;
}

