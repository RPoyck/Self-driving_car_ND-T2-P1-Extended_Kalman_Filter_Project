#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
#include "tools.h"

class FusionEKF {
    
    public:
	/**
	* Constructor.
	*/
	FusionEKF();

	/**
	* Destructor.
	*/
	virtual ~FusionEKF();

	/**
	* Run the whole flow of the Kalman Filter from here.
	*/
	void ProcessMeasurement(const MeasurementPackage &measurement_pack);

	/**
	* Kalman Filter update and prediction math lives in here.
	*/
	KalmanFilter ekf_;

    private:
	
	unsigned int min_measurements_for_outliers_;
	double maximum_include_distance_; 
	
	// check whether the tracking toolbox was initialised or not (first measurement)
	bool is_initialized_;
	unsigned long int no_of_measurements_;
	
	// previous time-stamp
	long long previous_timestamp_;

	// tool object used to compute Jacobian and RMSE
	Tools tools;
	Eigen::MatrixXd R_laser_;
	Eigen::MatrixXd R_radar_;
	Eigen::MatrixXd H_laser_;
	Eigen::MatrixXd Hj_;
    
	//acceleration noise components
	float noise_ax_;
	float noise_ay_;
	
	bool OutlierCheck(VectorXd x, const MeasurementPackage &measurement);

};

#endif /* FusionEKF_H_ */
