//-----------------------------------------------------------------------------
//
// This header file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Andrea Meo 2022. All rights reserved.
//
//-----------------------------------------------------------------------------

// Program headers
#ifndef HAMR_H_
#define HAMR_H_

// System headers
#include <string>
#include <vector>

namespace hamr{

   //-----------------------------------------------------------------------------
   // Function to initialise HAMR calculation
   //-----------------------------------------------------------------------------
   void initialize(const double Hmin,
                   const double Hmax,
                   const double system_dimensions_x,
                   const double system_dimensions_y,
                   const double system_dimensions_z,
                   const std::vector<double>& atom_coords_x,
                   const std::vector<double>& atom_coords_y,
                   const std::vector<double>& atom_coords_z,
                   const std::vector<int>& atom_type_array,
                   const int num_local_atoms);

   //-----------------------------------------------------------------------------
   // Function to check local temperature pulse is enabled
   //-----------------------------------------------------------------------------
	void fields(const int start_index,
					const int end_index,
					double H_applied,
					const double temperature,
					const double Tmin,
					const double Tmax,
					const double Hvecx,
					const double Hvecy,
					const double Hvecz,
					std::vector<double>& x_total_external_field_array,
					std::vector<double>& y_total_external_field_array,
					std::vector<double>& z_total_external_field_array);

   void hamr_continuous();
   // void hamr_pulsed();

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for hamr settings
   //-----------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //-----------------------------------------------------------------------------
   // Function to HAMR parameters
   //-----------------------------------------------------------------------------
   double get_head_position_x();
   double get_head_position_y();
   double get_head_speed();
   double get_FWHM_x();
   double get_FWHM_y();
   double get_laser_sigma_x();
   double get_laser_sigma_y();
   double get_field_bounds_x();
   double get_field_bounds_y();
   double get_bit_size();
   double get_track_size();
   double get_field_ramp_time();
   double get_NPS();
   bool get_initialisation_state();

   //-----------------------------------------------------------------------------------------------
   // Externally visible variables used for HAMR calculation
   //-----------------------------------------------------------------------------------------------
   extern bool head_laser_on;
   extern bool run_program;
   

} // end of hamr namespace

#endif // HAMR_H_
