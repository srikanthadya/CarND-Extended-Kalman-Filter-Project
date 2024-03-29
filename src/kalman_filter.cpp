#include "kalman_filter.h"
#define PI 3.14159265
using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_*x_;
  MatrixXd FT_ = F_.transpose();
  P_ = F_*P_*FT_ + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd y = z-H_*x_;
  MeasurementUpdate(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];
  
  float rho = sqrt(px*px + py*py);
  float theta = atan2(py, px);
  float rho_dot = (px*vx + py*vy)/rho;
  
  VectorXd h = VectorXd(3);
  h << rho,theta,rho_dot;
  VectorXd y = z-h;
  if( y[1] > PI )
    y[1] -= 2*PI;
  if( y[1] < -PI )
    y[1] += 2*PI;
  MeasurementUpdate(y);
  
}

void KalmanFilter::MeasurementUpdate(const VectorXd &y) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  MatrixXd HT_ = H_.transpose();
  MatrixXd S = H_*P_*HT_ + R_;
  MatrixXd Sinv = S.inverse();
  MatrixXd K_ = P_*HT_*Sinv;
  x_ = x_+ K_*y;
  int xsize = x_.size();
  MatrixXd I = MatrixXd::Identity(xsize,xsize);
  P_ = (I- K_*H_)*P_;
}