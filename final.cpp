#include <iostream>
#include <fstream>
#include <filesystem>
#include "csv.h"
#include "gradient.h"
#include "matrix.h"

//read csv into matrix
Matrix read_csv_to_matrix(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file");
    }

    Csv csv(file);
    std::vector<std::vector<double>> data;
    std::string line;

    while (csv.getline(line)) {
        std::vector<double> row;
        for (int i = 0; i < csv.getnfield(); ++i) {
            row.push_back(std::stod(csv.getfield(i)));
        }
        data.push_back(row);
    }

    int rows = data.size();
    int columns = data.empty() ? 0 : data[0].size();
    Matrix matrix(rows, columns);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            matrix.set_value(i, j, data[i][j]);
        }
    }

    return matrix;
}

int main() {
    std::string outputFolder = "output";
    //create the output directory
    //make sure /std:c++17 or higher is used when compiling
    std::filesystem::create_directory(outputFolder);

    try {
        Matrix returns = read_csv_to_matrix("asset_returns.csv");
        
        for (int i = 1; i < 21; ++i) {
            //increment the rate by 0.5%
            double returnRate = i * 0.005;
            std::string outputFilename = outputFolder + "/returns_" + std::to_string(returnRate) + ".csv";

            ConjugateGradientMethod conj_grad = ConjugateGradientMethod(returns);
            //train-test split is 100 days train, 12 day test
            Matrix result = conj_grad.train_test(100, 12, returnRate, 1e-6);
            result.to_csv(outputFilename, ",");
        }
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}