#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

// Constructor //
KalmanFilter::KalmanFilter() {}


// Destructor //
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
    TODO:
	* predict the state
    */
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}


void KalmanFilter::UpdateCore(const VectorXd &y) {
    
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
    
}


void KalmanFilter::Update(const VectorXd &z) {
    /**
    TODO:
    * update the state by using Kalman Filter equations
    */
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    
    UpdateCore(y);
    
}


VectorXd KalmanFilter::Normalise_y(const VectorXd &y) {
    VectorXd y_norm = y;

    while (y_norm(1) < -M_PI) {y_norm(1) += 2*M_PI;}
    while (M_PI < y_norm(1)) {y_norm(1) -= 2*M_PI;}
    
    return y_norm;
}


void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /**
    TODO:
	* update the state by using Extended Kalman Filter equations
    */
    
    // Define the equations which map the predicted location x' from Cartesian to polar coordinates //
    VectorXd h = VectorXd(3);
    h(0) = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
    h(1) = atan2(x_(1), x_(0));
    h(2) = (x_(0)*x_(2) + x_(1)*x_(3))/h(0);
    
    VectorXd y = z - h;
    
    VectorXd y_norm = Normalise_y(y);
    
    UpdateCore(y_norm);
    
}
