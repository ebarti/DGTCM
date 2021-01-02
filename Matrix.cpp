//
// Created by Eloi on 12/20/20.
//

#include "Matrix.h"
#include <math.h>

Matrix::Matrix(unsigned int xSize, unsigned int ySize, unsigned int zSize, double initialValue):
        m_size_x(xSize),
        m_size_y(ySize),
        m_size_z(zSize) {
    // Flat[x + WIDTH * (y + DEPTH * z)] = Original[x, y, z]
    m_matrix.resize(m_size_x * m_size_y * m_size_z, initialValue);
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

double Matrix::AbsMax() const {
    double max = 0.0;
    for (unsigned idx = 0; idx < m_size_x * m_size_y * m_size_z; idx++) {
        if(abs(m_matrix[idx]) > max) max = abs(m_matrix[idx]);
    }
    return max;
}

double Matrix::Max() const {
    double max = 0.0;
    for (unsigned idx = 0; idx < m_size_x * m_size_y * m_size_z; idx++) {
        if(m_matrix[idx] > max) max = abs(m_matrix[idx]);
    }
    return max;
}


Matrix Matrix::operator+(Matrix & A) {
    Matrix B(m_size_x, m_size_y, m_size_z, 0.0);
    for (unsigned idx = 0; idx < m_size_x * m_size_y * m_size_z; idx++) {
        B(idx) = m_matrix[idx] + A(idx);
    }
    return B;
}

Matrix Matrix::operator-(Matrix & A) {
    Matrix B(m_size_x, m_size_y, m_size_z, 0.0);
    for (unsigned idx = 0; idx < m_size_x * m_size_y * m_size_z; idx++) {
        B(idx) = m_matrix[idx] - A(idx);
    }
    return B;
}


Matrix Matrix::operator+(double A) {
    Matrix B(m_size_x, m_size_y, m_size_z, 0.0);
    for (unsigned idx = 0; idx < m_size_x * m_size_y * m_size_z; idx++) {
        B(idx) = m_matrix[idx] + A;
    }
    return B;
}

Matrix Matrix::operator-(double A) {
    Matrix B(m_size_x, m_size_y, m_size_z, 0.0);
    for (unsigned idx = 0; idx < m_size_x * m_size_y * m_size_z; idx++) {
        B(idx) = m_matrix[idx] - A;
    }
    return B;
}

Matrix Matrix::operator*(double A) {
    Matrix B(m_size_x, m_size_y, m_size_z, 0.0);
    for (unsigned idx = 0; idx < m_size_x * m_size_y * m_size_z; idx++) {
        B(idx) = m_matrix[idx] * A;
    }
    return B;
}

Matrix Matrix::operator/(double A) {
    Matrix B(m_size_x, m_size_y, m_size_z, 0.0);
    for (unsigned idx = 0; idx < m_size_x * m_size_y * m_size_z; idx++) {
        B(idx) = m_matrix[idx] / A;
    }
    return B;
}

double &Matrix::operator()(const unsigned int & x, const unsigned int & y, const unsigned int & z) {
    return m_matrix[x + m_size_x * (y + m_size_y * z)];
}

double Matrix::at(const unsigned int & x, const unsigned int & y, const unsigned int & z) const{
    return m_matrix[x + m_size_x * (y + m_size_y * z)];
}

double Matrix::at(const unsigned int & idx) const{
    return m_matrix[idx];
}


void Matrix::Copy(Matrix *otherMatrix) {
    for (unsigned idx = 0; idx < m_size_x * m_size_y * m_size_z; idx++) {
        this->m_matrix[idx] = otherMatrix->at(idx);
    }
}

double &Matrix::operator()(const unsigned int & idx) {
    return m_matrix[idx];
}



