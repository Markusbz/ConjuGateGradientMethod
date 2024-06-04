#ifndef MATRIX_H
#define MATRIX_H

class Matrix {
private:
    int rows;
    int columns;
    double **data;

public:
    // constructor for matrix creation
    Matrix(int rows_, int columns_);
    Matrix(const Matrix &other);
    Matrix& operator=(const Matrix &other);
    
    // destructor
    ~Matrix();

    //declare functions

    // displaying matrix; const member function
    void display() const;
    //get and set value
    double get_value(int row, int col);
    void set_value(int row, int col, double value);

    Matrix transpose() const;

    //matrix multiplication
    Matrix mult(const Matrix &other) const;
    //matrix inversion
    Matrix invert() const;
    //slicing the matrix
    Matrix slice(int row_start = 0, int row_end = -1, int col_start = 0, int col_end = -1) const;
    //sum function
    Matrix sum(int axis = 0) const;
    //product fucntion
    Matrix prod(int axis = 0) const;
    //mean function
    Matrix mean(int axis = 0) const;
    //create covariance matrix
    Matrix cov(int axis = 0, bool sample = true) const;

    //overloaded operators
    Matrix operator+(const Matrix &other) const;
    Matrix operator-(const Matrix &other) const;
    Matrix operator*(const Matrix &other) const;
    Matrix operator/(const Matrix &other) const;
    Matrix operator+(double value) const;
    Matrix operator-(double value) const;
    Matrix operator*(double value) const;
    Matrix operator/(double value) const;

    //dynamic row and column adjustments
    void add_row(const double *new_row);
    void add_column(const double *new_column);

    //get shape of the matrix
    std::pair<int, int> get_shape() const;

    //method to take the log of all values in matrix
    Matrix log() const;

    //transform single row/column matrices to an array
    double* to_array() const;

    //write amtrix to csv
    void to_csv(const std::string& filename, std::string delimiter) const;

};

#endif