#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  
  if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
    std::cout << "Inputs are not valid!" << std::endl;
    return rmse;
  }
      
  std::vector<VectorXd> residuals;
  for (size_t i = 0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    
    residual = residual.array() * residual.array();
    rmse += residual;
  }
      
  rmse = rmse / estimations.size();
  
  rmse = rmse.array().sqrt();
      
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  
  MatrixXd Hj(3, 4);
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  
  double denominator_square = px * px + py * py;
  double denominator = sqrt(denominator_square);
  double denominator_cube = denominator * denominator_square;
  
  if (denominator_square < 0.0000001) {
    std::cout << "Division by zero" << std::endl;
  } else {
    Hj(0, 0) = px / denominator;
    Hj(0, 1) = py / denominator;
    Hj(1, 0) = -py / denominator_square;
    Hj(1, 1) = px / denominator_square;
    Hj(2, 0) = py * (vx * py - vy * px) / denominator_cube;
    Hj(2, 1) = px * (vy * px - vx * py) / denominator_cube;
    Hj(2, 2) = px / denominator;
    Hj(2, 3) = py / denominator;
  }
  
  return Hj;
}
