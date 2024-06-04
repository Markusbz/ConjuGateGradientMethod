#include <iostream>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <fstream>
#include "matrix.h"

// constructor initializing rows and columns and allocating data
Matrix::Matrix(int rows_, int columns_) : rows(rows_), columns(columns_) {
    if (rows <= 0 || columns <= 0) {
        throw std::invalid_argument("Matrix dimensions must be greater than zero.");
    }
    data = new double*[rows]; // allocate memory for row pointers
    for (int i = 0; i < rows; i++) {
        data[i] = new double[columns]; // allocate memory for each row
        std::fill(data[i], data[i] + columns, 0.0); // Initialize all elements to 0
    }
}

// copy constructor for deep copying matrix data
Matrix::Matrix(const Matrix &other) : rows(other.rows), columns(other.columns) {
    data = new double*[rows];
    for (int i = 0; i < rows; i++) {
        data[i] = new double[columns];
        for (int j = 0; j < columns; j++) {
            data[i][j] = other.data[i][j]; // copy each element
        }
    }
}

// assignment operator to handle deep copy
Matrix& Matrix::operator=(const Matrix &other) {
    if (this == &other) {
        return *this; // handle self-assignment
    }

    // delete existing data
    for (int i = 0; i < rows; i++) {
        delete[] data[i];
    }
    delete[] data;

    // reallocate memory and copy data from the source
    rows = other.rows;
    columns = other.columns;
    data = new double*[rows];
    for (int i = 0; i < rows; i++) {
        data[i] = new double[columns];
        for (int j = 0; j < columns; j++) {
            data[i][j] = other.data[i][j];
        }
    }

    return *this;
}

//destructor
Matrix::~Matrix(){
    for(int i=0; i<rows;i++){
       delete[] data[i];
    }
    delete[] data;
}

//get a single value from the matrix
double Matrix::get_value(int row, int col) {
    if (row >= 0 && col>=0 && row<rows && col<columns){
        return data[row][col];
    }
    else{
        std::ostringstream oss;
        oss << "Index out of range: (" << row << ", " << col << ")";
        throw std::out_of_range(oss.str());
    }
}
//add a value to the matrix
void Matrix::set_value(int row, int col, double value){
    if (row >= 0 && col>=0 && row<rows && col<columns){
        data[row][col] = value;
    }
    else{
        std::ostringstream oss;
        oss << "Index out of range: (" << row << ", " << col << ")";
        std::cout << oss.str() << std::endl;
        throw std::out_of_range(oss.str());
    }
}
//display the matrix
void Matrix::display() const{
    for(int i=0;i<rows;i++){
        for(int j=0;j<columns;j++){
            std::cout << data[i][j] << " ";
        }
    std::cout << "" <<std::endl;
    }
}

// transpose matrix
Matrix Matrix::transpose() const {
    Matrix temp(columns, rows); //create matrix with inverted dimensions
    for(int row = 0; row < rows; row++){
        for(int col = 0; col < columns; col++){
            temp.set_value(col, row, data[row][col]);
        }
    }
    return temp;
}

Matrix Matrix::mult(const Matrix &other) const {
    if(columns==other.rows){
        //initialize new matrix
        Matrix temp(rows, other.columns);
        //iterate through columns and rows and update temporary sum
        for(int row=0; row < rows; row++) {
            for(int col=0; col < other.columns; col++){
                double temp_sum = 0;
                for(int i=0; i < columns; i++){
                    temp_sum += data[row][i]*other.data[i][col]; //compute dot product for each element
                }
                temp.data[row][col] = temp_sum;
            }
        }
        return temp;
    }
    else{
        std::ostringstream oss;
        oss << "Matrices with shape: (" << rows << ", " << columns << ") and (" << other.rows << ", " << other.columns << ") don't match";
        std::cout << oss.str() << std::endl;
        throw std::invalid_argument(oss.str());
    }
}

Matrix Matrix::invert() const {
    if (rows != columns) {
        throw std::invalid_argument("Matrix must be square to invert");
    }

    int n = rows;
    Matrix result(n, n);
    Matrix temp(*this);

    //initialize identity matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                result.set_value(i, j, 1.0);
            } else {
                result.set_value(i, j, 0.0);
            }
        }
    }

    //perform Gaussian elimination
    for (int i = 0; i < n; i++) {
        //find pivot
        double pivot = temp.get_value(i, i);
        if (pivot == 0) {
            throw std::runtime_error("Matrix is singular and cannot be inverted");
        }

        //scale pivot row
        for (int j = 0; j < n; j++) {
            temp.set_value(i, j, temp.get_value(i, j) / pivot);
            result.set_value(i, j, result.get_value(i, j) / pivot);
        }

        //eliminate other rows
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = temp.get_value(k, i);
                for (int j = 0; j < n; j++) {
                    temp.set_value(k, j, temp.get_value(k, j) - factor * temp.get_value(i, j));
                    result.set_value(k, j, result.get_value(k, j) - factor * result.get_value(i, j));
                }
            }
        }
    }

    return result;
}

Matrix Matrix::slice(int row_start, int row_end, int col_start, int col_end) const {
    if (row_end == -1) row_end = rows;
    if (col_end == -1) col_end = columns;

    if (row_start < 0 || row_start >= rows || row_end < 0 || row_end > rows || row_start >= row_end) {
        std::cout << "Row Index out of Range" << std::endl;
        throw std::out_of_range("Invalid slice range");
    }
    else {
        if(col_start < 0 || col_start >= columns || col_end < 0 || col_end > columns ||  col_start >= col_end) {
            std::cout << "Column Index out of Range" << std::endl;
            throw std::out_of_range("Invalid slice range");
        }

        else{
            int new_rows = row_end - row_start;
            int new_cols = col_end - col_start;
            Matrix result(new_rows, new_cols);

            for (int i = 0; i < new_rows; ++i) {
                for (int j = 0; j < new_cols; ++j) {
                    result.set_value(i, j, data[row_start + i][col_start + j]);
                }
            }
            return result;
        }
} }

//calculate the sum of elements along a specified axis
Matrix Matrix::sum(int axis) const{
    if(axis == 0){
        Matrix temp(1, columns); //sum across rows, result is 1 x columns
        for (int j = 0; j < columns; j++){
            double temp_sum = 0;
            for(int i=0; i< rows; i++){
                temp_sum += data[i][j];
            }
            temp.set_value(0, j, temp_sum);
        }
        return temp;
    }
    else if (axis == 1)
    {
        Matrix temp(rows, 1); //sum across rows, result is rows x 1
        for (int i = 0; i < rows; i++){
            double temp_sum = 0;
            for(int j=0;  j < columns; j++){
                temp_sum += data[i][j];
            }
            temp.set_value(i, 0, temp_sum);
        }
        return temp;
    }
    else{
        std::cout << "Key: Axis =" << axis << "not accepted" << std::endl;
        throw std::out_of_range("Axis key not accepted.");
    }
}

//calculate mean of elements along a specified axis
Matrix Matrix::mean(int axis) const{
    if(axis == 0){
        Matrix temp = sum(axis);
        for(int i=0; i < columns; i++){
            temp.set_value(0, i, temp.get_value(0, i)/rows);
        }
        return temp;
    }
    else if (axis == 1)
    {
        Matrix temp = sum(axis);
        for(int i=0; i < rows; i++){
            temp.set_value(i, 0, temp.get_value(i, 0)/columns);
        }
        return temp;
    }
    else{
        std::cout << "Axis key not accepted. Only 0 or 1." << std::endl;
        throw std::out_of_range("Axis key not accepted. Only 0 or 1.");
    }
}

// calculate the product of elements along a specified axis
Matrix Matrix::prod(int axis) const{
    if(axis == 0){
        Matrix temp(1, columns);
        for (int j = 0; j < columns; j++){
            double temp_sum = 1;
            for(int i=0; i< rows; i++){
                temp_sum *= data[i][j];
            }
            temp.set_value(0, j, temp_sum);
        }
        return temp;
    }
    else if (axis == 1)
    {
        Matrix temp(rows, 1);
        for (int i = 0; i < rows; i++){
            double temp_sum = 1;
            for(int j=0;  j < columns; j++){
                temp_sum *= data[i][j];
            }
            temp.set_value(i, 0, temp_sum);
        }
        return temp;
    }
    else{
        std::cout << "Key: Axis =" << axis << "not accepted" << std::endl;
        throw std::out_of_range("Axis key not accepted.");
    }
}

// calculate the covariance matrix along axis, sample or population
Matrix Matrix::cov(int axis, bool sample) const{
    if(axis==0){
        Matrix temp(columns, columns); //covariance matrix dimensions
        double divisor;
        if(sample){
            divisor = rows-1;
        }
        else{
            divisor = rows;
        }
        Matrix means = mean(0);
        //starting col
        for(int col=0; col < columns; col++){
            //covariance with col
            for(int col2 = col; col2 < columns; col2++){
                double temp_sum = 0;
                for(int i=0; i < rows; i++){
                    temp_sum += (data[i][col]-means.get_value(0, col))*(data[i][col2]-means.get_value(0, col2));
                }
                if(col != col2){
                    temp.set_value(col, col2, temp_sum/divisor); //fill up upper and lower "triangle" at the same time for faster performance
                    temp.set_value(col2, col, temp_sum/divisor);
                }
                else{
                    temp.set_value(col, col2, temp_sum/divisor);
                }
            }
        }
        return temp;
    }
    else if(axis==1){
        Matrix temp(rows, rows);
        double divisor;
        if(sample){
            divisor = columns-1;
        }
        else{
            divisor = columns;
        }
        Matrix means = mean(1);
        //starting row
        for(int row=0; row < rows; row++){
            //covariance with row
            for(int row2 = row; row2 < rows; row2++){
                double temp_sum = 0;
                for(int i=0; i < columns; i++){
                    temp_sum += (data[row][i]-means.get_value(row, 0))*(data[row2][i]-means.get_value(row2, 0));
                }
                if(row != row2){
                    temp.set_value(row, row2, temp_sum/divisor);
                    temp.set_value(row2, row, temp_sum/divisor);
                }
                else{
                    temp.set_value(row, row2, temp_sum/divisor);
                }
            }
        }
        return temp;
    }
    else{
        std::cout << "Key: Axis =" << axis << "not accepted" << std::endl;
        throw std::out_of_range("Axis key not accepted.");
    }
}


//overloaded operators
Matrix Matrix::operator+(const Matrix &other) const {
    if (rows == other.rows && columns == other.columns) {
        Matrix temp(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                temp.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return temp;
    } else {
        std::ostringstream oss;
        oss << "Matrices with shape: (" << rows << ", " << columns << ") and (" << other.rows << ", " << other.columns << ") don't match";
        throw std::invalid_argument(oss.str());
    }
}

Matrix Matrix::operator-(const Matrix &other) const {
    if (rows == other.rows && columns == other.columns) {
        Matrix temp(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                temp.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return temp;
    } else {
        std::ostringstream oss;
        oss << "Matrices with shape: (" << rows << ", " << columns << ") and (" << other.rows << ", " << other.columns << ") don't match";
        throw std::invalid_argument(oss.str());
    }
}

Matrix Matrix::operator*(const Matrix &other) const {
    if (rows == other.rows && columns == other.columns) {
        Matrix temp(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                temp.data[i][j] = data[i][j] * other.data[i][j];
            }
        }
        return temp;
    } else {
        std::ostringstream oss;
        oss << "Matrices with shape: (" << rows << ", " << columns << ") and (" << other.rows << ", " << other.columns << ") don't match";
        throw std::invalid_argument(oss.str());
    }
}

Matrix Matrix::operator/(const Matrix &other) const {
    if (rows == other.rows && columns == other.columns) {
        Matrix temp(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                if (other.data[i][j] == 0) {
                    throw std::invalid_argument("Division by zero in matrix element");
                }
                temp.data[i][j] = data[i][j] / other.data[i][j];
            }
        }
        return temp;
    } else {
        std::ostringstream oss;
        oss << "Matrices with shape: (" << rows << ", " << columns << ") and (" << other.rows << ", " << other.columns << ") don't match";
        throw std::invalid_argument(oss.str());
    }
}

Matrix Matrix::operator+(double value) const {
    Matrix temp(rows, columns);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            temp.data[i][j] = data[i][j] + value;
        }
    }
    return temp;
}

Matrix Matrix::operator-(double value) const {
    Matrix temp(rows, columns);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            temp.data[i][j] = data[i][j] - value;
        }
    }
    return temp;
}

Matrix Matrix::operator*(double value) const {
    Matrix temp(rows, columns);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            temp.data[i][j] = data[i][j] * value;
        }
    }
    return temp;
}

Matrix Matrix::operator/(double value) const {
    if (value == 0) {
        throw std::invalid_argument("Division by zero");
    }
    Matrix temp(rows, columns);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            temp.data[i][j] = data[i][j] / value;
        }
    }
    return temp;
}

void Matrix::add_row(const double *new_row) {
    // Create a new data array with an extra row
    double **new_data = new double*[rows + 1];
    for (int i = 0; i < rows + 1; i++) {
        new_data[i] = new double[columns];
    }

    // Copy the existing data to the new data array
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            new_data[i][j] = data[i][j];
        }
    }

    // Add the new row
    for (int j = 0; j < columns; j++) {
        new_data[rows][j] = new_row[j];
    }

    // Delete the old data array
    for (int i = 0; i < rows; i++) {
        delete[] data[i];
    }
    delete[] data;

    // Point data to the new data array and update rows
    data = new_data;
    rows++;
}

void Matrix::add_column(const double *new_column) {
    // Create a new data array with an extra column
    double **new_data = new double*[rows];
    for (int i = 0; i < rows; i++) {
        new_data[i] = new double[columns + 1];
    }

    // Copy the existing data to the new data array
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            new_data[i][j] = data[i][j];
        }
    }

    // Add the new column
    for (int i = 0; i < rows; i++) {
        new_data[i][columns] = new_column[i];
    }

    // Delete the old data array
    for (int i = 0; i < rows; i++) {
        delete[] data[i];
    }
    delete[] data;

    // Point data to the new data array and update columns
    data = new_data;
    columns++;
}

std::pair<int, int> Matrix::get_shape() const {
    return {rows, columns};
}

//logging all values in the matrix
Matrix Matrix::log() const {
    Matrix result(rows, columns);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            if (data[i][j] <= 0) {
                throw std::domain_error("Logarithm of non-positive value is undefined");
            }
            result.set_value(i, j, std::log(data[i][j]));
        }
    }
    return result;
}

//convert matrix to array
double* Matrix::to_array() const {
    if (rows != 1 && columns != 1) {
        throw std::invalid_argument("Matrix must be a single row or single column to convert to array");
    }
    double* array = new double[std::max(rows, columns)];
    if (rows == 1) {
        for (int i = 0; i < columns; ++i) {
            array[i] = data[0][i];
        }
    } else {
        for (int i = 0; i < rows; ++i) {
            array[i] = data[i][0];
        }
    }
    return array;
}

void Matrix::to_csv(const std::string& filename, std::string delimiter) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for writing");
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            file << data[i][j];
            if (j < columns - 1) file << delimiter;  //separate columns with comma
        }
        file << "\n";  // End of row
    }
    file.close();
}