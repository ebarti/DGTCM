//
// Created by Eloi on 12/20/20.
//

#include "Matrix.h"
#include <math.h>

Matrix::Matrix(unsigned int xSize, unsigned int ySize, unsigned int zSize, double initialValue):
        m_size_x(xSize),
        m_size_y(ySize),
        m_size_z(zSize) {
    m_matrix.resize(m_size_x);
    for (unsigned x = 0; x < m_matrix.size(); x++) {
        m_matrix[x].resize(m_size_y);
        for (unsigned y = 0; y < m_matrix[x].size(); y++) {
            m_matrix[x][y].resize(m_size_z, initialValue);
        }
    }
}

Matrix::Matrix(const Matrix & A):
        m_matrix(A.m_matrix),
        m_size_x(A.m_size_x),
        m_size_y(A.m_size_y),
        m_size_z(A.m_size_z){

}


void Matrix::print() const {

}

unsigned int Matrix::XSize() const {
    return m_size_x;
}

unsigned int Matrix::YSize() const {
    return m_size_y;
}

unsigned int Matrix::ZSize() const {
    return m_size_z;
}

double Matrix::Max() const {
    double max = 0.0;
    for (unsigned x = 0; x < m_size_x; x++) {
        for (unsigned y = 0; y < m_size_y; y++) {
            for (unsigned z = 0; z < m_size_z; z++) {
                if(abs(m_matrix[x][y][z]) > max) max = abs(m_matrix[x][y][z]);
            }
        }
    }
    return max;
}

double Matrix::AbsMax() const {
    double max = 0.0;
    for (unsigned x = 0; x < m_size_x; x++) {
        for (unsigned y = 0; y < m_size_y; y++) {
            for (unsigned z = 0; z < m_size_z; z++) {
                if(m_matrix[x][y][z] > max) max = m_matrix[x][y][z];
            }
        }
    }
    return max;
}


Matrix Matrix::operator+(Matrix & A) {
    Matrix B(m_size_x, m_size_y, m_size_z, 0.0);
    for (unsigned x = 0; x < m_size_x; x++) {
        for (unsigned y = 0; y < m_size_y; y++) {
            for (unsigned z = 0; z < m_size_z; z++) {
                B(x, y, z) = m_matrix[x][y][z] + A(x, y, z);
            }
        }
    }
    return B;
}

Matrix Matrix::operator-(Matrix & A) {
    Matrix B(m_size_x, m_size_y, m_size_z, 0.0);
    for (unsigned x = 0; x < m_size_x; x++) {
        for (unsigned y = 0; y < m_size_y; y++) {
            for (unsigned z = 0; z < m_size_z; z++) {
                B(x, y, z) = m_matrix[x][y][z] - A(x, y, z);
            }
        }
    }
    return B;
}


Matrix Matrix::operator+(double A) {
    Matrix B(m_size_x, m_size_y, m_size_z, 0.0);
    for (unsigned x = 0; x < m_size_x; x++) {
        for (unsigned y = 0; y < m_size_y; y++) {
            for (unsigned z = 0; z < m_size_z; z++) {
                B(x, y, z) = m_matrix[x][y][z] + A;
            }
        }
    }
    return B;
}

Matrix Matrix::operator-(double A) {
    Matrix B(m_size_x, m_size_y, m_size_z, 0.0);
    for (unsigned x = 0; x < m_size_x; x++) {
        for (unsigned y = 0; y < m_size_y; y++) {
            for (unsigned z = 0; z < m_size_z; z++) {
                B(x, y, z) = m_matrix[x][y][z] - A;
            }
        }
    }
    return B;
}

Matrix Matrix::operator*(double A) {
    Matrix B(m_size_x, m_size_y, m_size_z, 0.0);
    for (unsigned x = 0; x < m_size_x; x++) {
        for (unsigned y = 0; y < m_size_y; y++) {
            for (unsigned z = 0; z < m_size_z; z++) {
                B(x, y, z) = m_matrix[x][y][z] * A;
            }
        }
    }
    return B;
}

Matrix Matrix::operator/(double A) {
    Matrix B(m_size_x, m_size_y, m_size_z, 0.0);
    for (unsigned x = 0; x < m_size_x; x++) {
        for (unsigned y = 0; y < m_size_y; y++) {
            for (unsigned z = 0; z < m_size_z; z++) {
                B(x, y, z) = m_matrix[x][y][z] / A;
            }
        }
    }
    return B;
}

double &Matrix::operator()(const unsigned int & x, const unsigned int & y, const unsigned int & z) {
    return m_matrix[x][y][z];
}

double Matrix::at(const unsigned int & x, const unsigned int & y, const unsigned int & z) const{
    return m_matrix[x][y][z];;
}

double Matrix::MaxDiff(const Matrix *&iOther) const {
    double maxDiff = 0.0;
    for(unsigned int x = 0; x < m_size_x; x++) {
        for(unsigned int y = 0; y < m_size_y; y++) {
            for(unsigned int z = 0; z < m_size_z; z++) {
                double val = abs(this->at(x,y,z) - iOther->at(x,y,z));
                if(val > maxDiff) maxDiff = val;
            }
        }   
    }
    return maxDiff;
}

void Matrix::Copy(Matrix *otherMatrix) {
    for(unsigned int x = 0; x < m_size_x; x++) {
        for(unsigned int y = 0; y < m_size_y; y++) {
            for(unsigned int z = 0; z < m_size_z; z++) {
                this->m_matrix[x][y][z] = otherMatrix->at(x,y,z);
            }
        }
    }
}



