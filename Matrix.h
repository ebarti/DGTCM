//

// Created by Eloi on 12/20/20.
//

#ifndef DGTCM_MATRIX_H
#define DGTCM_MATRIX_H
#include <vector>
#include <ostream>

class Matrix {
public:
    Matrix(unsigned int xSize, unsigned int ySize, unsigned int zSize, double initialValue);
    Matrix(const Matrix&);

    // Public Methods
    void print() const;
    unsigned int XSize() const;
    unsigned int YSize() const;
    unsigned int ZSize() const;
    double Max() const;
    double AbsMax() const;
    void Copy(Matrix * otherMatrix);
    // Operators
    Matrix operator+(Matrix &);
    Matrix operator-(Matrix &);
    Matrix operator+(double);
    Matrix operator-(double);
    Matrix operator*(double);
    Matrix operator/(double);
    double& operator()(const unsigned &, const unsigned &, const unsigned &);
    double& operator()(const unsigned &);
    double at(const unsigned &, const unsigned &, const unsigned &) const;
private:
    unsigned m_size_x;
    unsigned m_size_y;
    unsigned m_size_z;

    std::vector<double> m_matrix;
};


#endif //DGTCM_MATRIX_H
