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
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

    
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
        float rho = measurement_pack.raw_measurements_[0];
        float phi = measurement_pack.raw_measurements_[1];
        float rho_dot = measurement_pack.raw_measurements_[2];
 
        float x = rho * cos(phi);
        float y = rho * sin(phi);
        float vx = rho_dot * cos(phi);
        float vy = rho_dot * sin(phi);
        ekf_.x_ << x, y, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
        ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
      float verySmallNumber = 0.0001;
      if (fabs(ekf_.x_(0)) < verySmallNumber and fabs(ekf_.x_(1)) < verySmallNumber){
          ekf_.x_(0) = verySmallNumber;
          ekf_.x_(1) = verySmallNumber;
      }
      
      ekf_.P_ = MatrixXd(4, 4);
      ekf_.P_ << 1, 0, 0, 0,
			   0, 1, 0, 0,
			   0, 0, 200, 0,
			   0, 0, 0, 200;
      previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

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
    float dt = (measurement_pack.timestamp_ - previous_timestamp_);
    dt /= 1000000.0; // convert micros to s
    previous_timestamp_ = measurement_pack.timestamp_;
    
    // State transition matrix update
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, dt, 0,
               0, 1, 0, dt,
               0, 0, 1, 0,
               0, 0, 0, 1;
    // Noise covariance matrix computation
    // Noise values from the task
    float noise_ax = 9.0;
    float noise_ay = 9.0;
    // Precompute some usefull values to speed up calculations of Q
    float dt2 = dt * dt; //dt^2
    float dt3 = dt2 * dt; //dt^3
    float dt4 = dt3 * dt; //dt^4
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt4/4 * noise_ax, 0, dt3/2 * noise_ax, 0,
               0, dt4/4 * noise_ay, 0, dt3/2 * noise_ay,
               dt3/2 * noise_ax, 0, dt2 * noise_ax, 0,
               0, dt3/2 * noise_ay, 0, dt2 * noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
      
  } else {
    // Laser updates
      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
