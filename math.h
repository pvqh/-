#pragma once

#include <string>
namespace Matrix
{
    double** Add(double** matrix1, double** matrix2, int row, int col)
    {
        //double result[row][col]{};
        double** result = new double*[row];
        for (int i = 0; i < row; i++) {
            result[i] = new double[col];
        }
        for(int i = 0; i < row; i++)
        {
            for(int j = 0; j < col; j++)
            {
                result[i][j] = matrix1[i][j] + matrix2[i][j]; 
            }
        }
        return result;
    }
    double** Subtract(double** matrix1, double** matrix2, int row, int col)
    {
        double** result = new double*[row];
        for (int i = 0; i < row; i++) {
            result[i] = new double[col];
        }
        for(int i = 0; i < row; i++)
        {
            for(int j = 0; j < col; j++)
            {
                result[i][j] = matrix1[i][j] - matrix2[i][j]; 
            }
        }
        return result;
    }
    double** Multiply(double** matrix1, double** matrix2, int row1, int col1, int row2, int col2)
    {
        double** result = new double*[row1];
        for (int i = 0; i < row1; i++) {
            result[i] = new double[col2];
        }
        for(int i = 0; i < row1; i++)
        {
            for(int j = 0; j < col2; j++)
            {
                for(int k = 0; k < col1; k++)
                result[i][j] += matrix1[i][k] * matrix2[k][j]; 
            }
        }
        return result;
    }
    double Determinant(double** matrix, int size)
    {

    int det = 0;
    int sign = 1;
 
    // Base Case
    if (size == 1) {
        det = matrix[0][0];
    }
    else if (size == 2) {
        det = (matrix[0][0] * matrix[1][1])
              - (matrix[0][1] * matrix[1][0]);
    }
 
    // Perform the Laplace Expansion
    else {
        for (int i = 0; i < size; i++) {
 
            // Stores the cofactor matrix
            double** cofactor = new int*[size - 1];
            for (int j = 0; j < size - 1; j++) {
                cofactor[j] = new double[size - 1];
            }
            int sub_i = 0, sub_j = 0;
            for (int j = 1; j < size; j++) {
                for (int k = 0; k < size; k++) {
                    if (k == i) {
                        continue;
                    }
                    cofactor[sub_i][sub_j] = matrix[j][k];
                    sub_j++;
                }
                sub_i++;
                sub_j = 0;
            }
 
            // Update the determinant value
            det += sign * matrix[0][i]
                   * Determinant(cofactor, size - 1);
            sign = -sign;
            for (int j = 0; j < size - 1; j++) {
                delete[] cofactor[j];
            }
            delete[] cofactor;
        }
    }
 
    // Return the final determinant value
    return det;
    }
}

namespace Function
{
    struct Func
    {
        std::string f;
        double x1, x2, x3 = 0.0;
        double df_dx1 = -4*x1 + 3*x2 - x3 + 8;
        double df_dx2 = 3*x1 - 8*x2;
        double df_dx3 = -12*x3 - x1;

        int df_dx1x1 = -4;
        int df_dx1x2 = 3;
        int df_dx1x3 = -1;
        int df_dx2x1 = df_dx1x2;
        int df_dx2x2 = -8;
        int df_dx2x3 = 0;
        int df_dx3x1 = df_dx1x3;
        int df_dx3x2 = df_dx2x3;
        int df_dx3x3 = -12;

        double gradiant_f [3][1] = {
            {df_dx1}, {df_dx2}, {df_dx3}
        };

        int Jacobian[3][3] = {
            {df_dx1x1, df_dx1x2, df_dx1x3},
            {df_dx2x1, df_dx2x2, df_dx2x3},
            {df_dx3x1, df_dx3x2, df_dx3x3}
        };

    };
    double Partial_Derivative(double x0, double y0, double z0);
}