#include <stdexcept>
#include <cmath>
#include <iostream>
#include <sstream>
#include "gradient.h"
#include "matrix.h"

ConjugateGradientMethod::ConjugateGradientMethod(const Matrix &returns) : returns(returns), n(returns.get_shape().second) {}

Matrix ConjugateGradientMethod::build_Q() const{
    //calculate the covariance (sample cov) matrix (Σ) along columns
    Matrix Q = returns.cov(0, true);

    //calculate the mean returns (r)
    Matrix meanReturns = returns.mean(0)*-1;

    //create a vector of ones (e)
    Matrix onesVector(n, 1);
    for (int i = 0; i < n; ++i) {
        onesVector.set_value(i, 0, -1.0);
    }

    //add -r and -e as columns to Q
    Q.add_column(meanReturns.to_array());
    Q.add_column(onesVector.to_array());
    //add -r and -e as rows to Q
    Q.add_row(meanReturns.to_array());
    Q.add_row(onesVector.to_array());

    //set the last elements to 0
    Q.set_value(n, n, 0.0);
    Q.set_value(n, n + 1, 0.0);
    Q.set_value(n + 1, n, 0.0);
    Q.set_value(n + 1, n + 1, 0.0);

    return Q;
}

Matrix ConjugateGradientMethod::build_b(double target_return) {
        Matrix b(n+2,1);
        for(int i=0; i < n; i++){
            b.set_value(i, 0, 0.);
        }
        b.set_value(n, 0, target_return);
        b.set_value(n+1, 0, -1.);  //set the sum of weights to be 1

        return b;
    }

//initial guess for the optimization
Matrix ConjugateGradientMethod::build_x0() {
    double init_weight = 1./n;
    Matrix x0(n+2, 1);
    //set initial weightings for all assets equal
    for(int i = 0; i < n; i++){
        x0.set_value(i, 0, init_weight);
    }
    //set initial lagrange multipliers to 0
    x0.set_value(n, 0, 0.0);
    x0.set_value(n+1, 0, 0.0);

    return x0;
}

//run conjugate gradient method to find the optimal portfolio weights
Matrix ConjugateGradientMethod::run(double target_return, double tolerance) {
    Matrix Q = build_Q();
    Matrix x = build_x0();
    //target-return negative, as we want to minimize
    Matrix s = build_b(-target_return) - Q.mult(x);
    Matrix p = s;
    int iteration = 0;

    Matrix s1(n+2,1);
    double a;
    double beta;


    double e = (s.transpose().mult(s).get_value(0,0));
    while (e > tolerance) //iterate until e is within the tolerance
    {
        if (iteration >= 1000) {
            throw std::runtime_error("Conjugate gradient method failed to converge within " + std::to_string(1000) + " iterations.");
        }
        a = e/(p.transpose().mult(Q).mult(p).get_value(0,0));
        x = x + p*a;  //update solution x
        s = s - (Q*a).mult(p);
        beta = (s.transpose().mult(s).get_value(0, 0))/e;
        p = s + p*beta;
        std::cout << e << std::endl;

        e = (s.transpose().mult(s).get_value(0,0));

        iteration++;
    }
    return x.slice(0, n,0,-1); //return the first n elements of x, excluding lagrange multipliers
}

//train and test the model using rolling windows for parameter estimation and backtesting
Matrix ConjugateGradientMethod::train_test(int train_window, int test_window, double target_return, double tolerance){
    int rows = returns.get_shape().first; //total number of days for which data is available
    int size = (rows-train_window)/test_window; //calculate the number of testing periods

    Matrix result(size, n+2); //prepare a result matrix for weights, return, variance
    Matrix covar_matrix(n, n);

    double return_i;
    double portfolio_var;

    int i = train_window;
    int j;

    while((i+test_window) < rows){
        Matrix train_data = returns.slice((i-train_window), i, 0, -1);
        Matrix test_data = returns.slice(i, i+test_window, 0, -1);


        ConjugateGradientMethod train_method(train_data); //initialize the conjugate gradient method with training data
        Matrix weights = train_method.run(target_return, tolerance); //run optimization to get portfolio weights

        
        return_i = test_data.mult(weights).mean(0).get_value(0,0); //calculate average return using the trained weights
        covar_matrix = test_data.cov(0,true); //calculate covariance matrix for the test data

        portfolio_var = weights.transpose().mult(covar_matrix).mult(weights).get_value(0,0); //calculate portfolio variance

        j = (i-train_window)/test_window;

        for(int b = 0; b < n; b++){
            result.set_value(j, b, weights.get_value(b, 0));
        }
        result.set_value(j, n, return_i);
        result.set_value(j, n+1, portfolio_var);

        i += test_window;
    }
    return result; //return the matrix containing weights, returns, and variances for each period
}
