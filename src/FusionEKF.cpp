#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    
    min_measurements_for_outliers_ = 10;	// Minimum number of measurements before measurement outliers are discarded from the update step //
    maximum_include_distance_ = 2.0;		// [m] Maximum absolute distance from the current prediction at which measurements are included for the update step //
    
    is_initialized_ = false;
    no_of_measurements_ = 0;

    previous_timestamp_ = 0;

    // initialising matrices
    ekf_.x_ = VectorXd(4);

    /**
    TODO:
	* Finish initialising the FusionEKF.
	* Set the process and measurement noises
    */

    // State covariance matrix P //
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ <<	1, 	0, 	0, 	0,
		0, 	1, 	0, 	0,
		0, 	0, 	1000, 	0,
		0, 	0, 	0, 	1000;

    //measurement covariance matrix - laser
    R_laser_ = MatrixXd(2, 2);
    R_laser_ <<	0.0225, 0,
		0, 	0.0225;

    //measurement covariance matrix - radar
    R_radar_ = MatrixXd(3, 3);
    R_radar_ << 0.09, 	0, 	0,
		0, 	0.0009, 0,
		0, 	0, 	0.09;

    // Measurement matrix //
    H_laser_ = MatrixXd(2, 4);
    H_laser_ <<	1, 0, 0, 0,
		0, 1, 0, 0;
		
    // Initial H Jacobian matrix //
    Hj_ = MatrixXd(3, 4);
    Hj_ <<  1, 1, 0, 0,
	    1, 1, 0, 0,
	    1, 1, 1, 1;
		
    // The initial transition matrix F_ //
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 	1, 0, 1, 0,
		0, 1, 0, 1,
		0, 0, 1, 0,
		0, 0, 0, 1;
	    
    //set the acceleration noise components

    noise_ax_ = 9;
    noise_ay_ = 9;
	
}


/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}


bool FusionEKF::OutlierCheck(VectorXd x, const MeasurementPackage &measurement_pack) {
    
    bool is_outlier = false;
    
    if (no_of_measurements_ < min_measurements_for_outliers_) {
	no_of_measurements_++;
	return is_outlier;
    }
    
    double x_value, y_value;
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	double rho = measurement_pack.raw_measurements_[0];
	double phi = measurement_pack.raw_measurements_[1];
	x_value = rho * cos(phi);
	y_value = rho * sin(phi);
    }
    else {
	if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
	    x_value = measurement_pack.raw_measurements_[0];
	    y_value = measurement_pack.raw_measurements_[1];
	}
    }
    
    double d_x = abs(x(0) - x_value);
    double d_y = abs(x(1) - y_value);
    
    if (maximum_include_distance_ < d_x || maximum_include_distance_ < d_y ) {
	cout << "It's an outlier! dx: " << d_x << "  dy: " << d_y << "  Type: " << measurement_pack.sensor_type_ << endl;
	no_of_measurements_--;
	is_outlier = true;
    }
    else {no_of_measurements_++;}
    
    return is_outlier;
    
}


void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    
    /*****************************************************************************
    *  Initialisation
    ****************************************************************************/
    if (!is_initialized_) {
	/**
	TODO:
	* Initialise the state ekf_.x_ with the first measurement.
	* Create the covariance matrix.
	* Remember: you'll need to convert radar from polar to Cartesian coordinates.
	*/
	// first measurement
	cout << "EKF: " << endl;
	ekf_.x_ = VectorXd(4);

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	    /**
	    Convert radar from polar to Cartesian coordinates and initialise state.
	    */
	    double rho = measurement_pack.raw_measurements_[0];
	    double phi = measurement_pack.raw_measurements_[1];
	    double rhoDot = measurement_pack.raw_measurements_[2];
	    
	    ekf_.x_(0) = rho * cos(phi);
	    ekf_.x_(1) = rho * sin(phi);
	    ekf_.x_(2) = rhoDot * cos(phi);
	    ekf_.x_(3) = rhoDot * sin(phi);
	}
	else {
	    if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
		/**
		Initialise state.
		*/
		//set the state with the initial location and zero velocity
		ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
	    }
	}

	previous_timestamp_ = measurement_pack.timestamp_;
	
	// done initialising, no need to predict or update
	is_initialized_ = true;
	return;
    }

    //compute the time elapsed between the current and previous measurements
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds

    previous_timestamp_ = measurement_pack.timestamp_;
    
    /*****************************************************************************
    *  Prediction
    ****************************************************************************/
	
    /**
    TODO:
	* Update the state transition matrix F according to the new elapsed time.
	- Time is measured in seconds.
	* Update the process noise covariance matrix.
	* Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
    */
    
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;

    //Modify the F matrix so that the time is integrated

    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;
    
    //set the process covariance matrix Q

    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ <<  dt_4/4*noise_ax_, 	0, 			dt_3/2*noise_ax_, 	0,
		0, 			dt_4/4*noise_ay_,	0, 			dt_3/2*noise_ay_,
		dt_3/2*noise_ax_, 	0, 			dt_2*noise_ax_,		0,
		0, 			dt_3/2*noise_ay_,	0, 			dt_2*noise_ay_;
    
    ekf_.Predict();

    /*****************************************************************************
    *  Update
    ****************************************************************************/

    /**
    TODO:
	* Use the sensor type to perform the update step.
	* Update the state and covariance matrices.
    */
    
    if (!OutlierCheck(ekf_.x_, measurement_pack))
    {
	// When updating with Radar //
	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	    Tools tools;
	    ekf_.R_ = R_radar_;
	    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	} 
	// When updating with Lidar //
	else {
	    // Laser updates
	    ekf_.R_ = R_laser_;
	    ekf_.H_ = H_laser_;
	    ekf_.Update(measurement_pack.raw_measurements_);
	}
    }

    // print the output
//     cout << "x_ = " << ekf_.x_ << endl;
//     cout << "P_ = " << ekf_.P_ << endl;
}
