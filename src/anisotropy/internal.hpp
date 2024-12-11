//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Jack Collings, Roberto Moreno Ortega, Sam Westmoreland and
//   Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef ANISOTROPY_INTERNAL_H_
#define ANISOTROPY_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the anisotropy module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <vector>

// Vampire headers
#include "material.hpp"
#include "anisotropy.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

   namespace internal{

      //-----------------------------------------------------------------------------
      // Internal data type definitions
      //-----------------------------------------------------------------------------

      //--------------------------------------------------------------------
      // Class to contain parameters for lattice anisotropy calculation.
      //
      // Tabulated values are read from a file and added point-wise to
      // the class. During variable initialisation interpolating functions
      // are determined to calculate k(T)
      //
      class lattice_anis_t
      {

         private:

            unsigned int Tmax; // maximum array value in tabulated function
            double k_Tmax; // value of anisotropy at k_Tmax (used for all T > Tmax)

            std::vector< unsigned int > T; // input temperature points
            std::vector< double > k; // input lattice anisotropy points
            std::vector< double > m; // calculated m value
            std::vector< double > c; // calculated c value

         public:

            void add_point( double, double );
            void set_interpolation_table();
            double get_lattice_anisotropy_constant( double );
            void output_interpolated_function( int );

      };

      //---------------------------------
      // struct for storing unit vectors
      //---------------------------------
      struct evec_t
      {
         double x;
         double y;
         double z;
      };

      //-----------------------------------------------------------------------------
      // materials class for storing anisotropy material parameters
      //-----------------------------------------------------------------------------
      class mp_t
      {

         private:

         public:

            // variables
            double ku2;       // second order uniaxial anisotropy constant ( - Ku1 )
            double k2r1;      // second order theta first order phi anisotropy constant
            double k2r1_odd;  // second order theta first order phi odd anisotropy constant
            double k2r2;      // second order theta second order phi anisotropy constant
            double k2r2_odd;  // seconf order theta second order phi odd anisotropy constant
            double ku4;       // fourth order uniaxial anisotropy constant ( ~Ku2 )
            double k4r1;      // fourth order theta first order phi anisotropy constant
            double k4r1_odd;  // fourth order theta first order phi odd anisotropy constant
            double k4r2;      // fourth order theta second order phi anisotropy constant
            double k4r2_odd;  // fourth order theta second order phi odd anisotropy constant
            double k4r3;      // fourth order theta third order phi odd anisotropy constant
            double k4r3_odd;  // fourth order theta third order phi odd anisotropy constant
            double k4r4;      // fourth order theta fourth order phi anisotropy constant
            double k4r4_odd;  // fourth order theta fourth order phi odd anisotropy constant
            double ku6;       // sixth order uniaxial anisotropy constant  ( ~Ku3 )
            double k6r1;      // sixth order theta first order phi anisotropy constant
            double k6r1_odd;  // sixth order theta first order phi odd anisotropy constant
            double k6r2;      // sixth order theta second order phi anisotropy constant
            double k6r2_odd;  // sixth order theta second order phi odd anisotropy constant
            double k6r3;      // sixth order theta third order phi anisotropy constant
            double k6r3_odd;  // sixth order theta third order phi odd anisotropy constant
            double k6r4;      // sixth order theta fourth order phi anisotropy constant
            double k6r4_odd;  // sixth order theta fourth order phi odd anisotropy constant
            double k6r5;      // sixth order theta fifth order phi anisotropy constant
            double k6r5_odd;  // sixth order theta fifth order phi odd anisotropy constant
            double k6r6;      // sixth order theta sixth order phi anisotropy constant
            double k6r6_odd;  // sixth order theta sixth order phi odd anisotropy constant

            double kc4; // fourth order cubic anisotropy constant (Kc1)
            double kc6; // sixth order cubic anisotropy constant (Kc2)

            double k_lattice; // uniaxial lattice anisotropy constant

            std::vector< double > kij; // surface/Neel anisotropy pair constant

            std::vector< double > ku_vector; // unit vector defining axis for uniaxial anisotropy
            std::vector< double > kr_vector; // unit vector defining axis for rotational anisotropy
            std::vector< double > kl_vector; // last unit vector perpendicular to kr_vector and ku_vector

            std::vector< double > u1_vector; // unit vector defining axis for uniaxial anisotropy
            std::vector< double > u2_vector; // unit vector defining axis for uniaxial anisotropy                        

            std::vector< double > kc_vector1; // first unit vector defining axis for cubic anisotropy
            std::vector< double > kc_vector2; // second unit vector defining axis for cubic anisotropy
            std::vector< double > kc_vector3; // third unit vector defining axis for cubic anisotropy

            std::vector< double > ku_tensor; // uniaxial second order anisotropy tensor
            std::vector< double > kc_tensor; // cubic fourth order anisotropy tensor

            lattice_anis_t lattice_anisotropy; // class containing lattice anisotropy data

            bool random_anisotropy; // flag to control random anisotropy by material
      		bool random_grain_anisotropy; // flag to control random anisotropy by grain

            // constructor
            mp_t ( const unsigned int max_materials = 666 ):
            	ku2( 0.0 ),       // set initial value of ku2 to zero
               k2r1( 0.0 ),      // set initial value of k2r1 to zero
               k2r1_odd( 0.0 ),  // set initial value of k2r1_odd to zero
               k2r2( 0.0 ),      // set initial value of k2r2 to zero
               k2r2_odd( 0.0 ),  // set initial value of k2r2_odd to zero
               ku4( 0.0 ),       // set initial value of ku4 to zero
               k4r1( 0.0 ),      // set initial value of k4r1 to zero
               k4r1_odd( 0.0 ),  // set initial value of k4r1_odd to zero
               k4r2( 0.0 ),      // set intital value of k4r2 to zero
               k4r2_odd( 0.0 ),  // set initial value of k4r2_odd to zero
               k4r3( 0.0 ),      // set initial value of k4r3_odd to zero
               k4r3_odd( 0.0 ),  // set initial value of k4r3_odd to zero
               k4r4( 0.0 ),      // set initial value of k4r4 to zero
               k4r4_odd( 0.0 ),  // set initial value of k4r4_odd to zero
               ku6( 0.0 ),       // set initial value of ku6 to zero
               k6r1( 0.0 ),      // set initial value of k6r1 to zero
               k6r1_odd( 0.0 ),  // set initial value of k6r1_odd to zero
               k6r2( 0.0 ),      // set initial value of k6r2 to zero
               k6r2_odd( 0.0 ),  // set initial value of k6r2_odd to zero
               k6r3( 0.0 ),      // set initial value of k6r3 to zero
               k6r3_odd( 0.0 ),  // set initial value of k6r3_odd to zero
               k6r4( 0.0 ),      // set initial value of k6r4 to zero
               k6r4_odd( 0.0 ),  // set initial value of k6r4_odd to zero
               k6r5( 0.0 ),      // set initial value of k6r5 to zero
               k6r5_odd( 0.0 ),  // set initial value of k6r5_odd to zero
               k6r6( 0.0 ),      // set initial value of k6r6 to zero
               k6r6_odd( 0.0 ),  // set initial value of k6r6_odd to zero

               kc4( 0.0 ), // set initial value of kc4 to zero
               kc6( 0.0 ), // set initial value of kc6 to zero

               k_lattice( 0.0 ), // set initial value of k_lattice to zero
               random_anisotropy( false ), // disable random anisotropy
               random_grain_anisotropy( false ) // disable random grain anisotropy
            {
               // resize arrays to correct size
               kij.resize( max_materials, 0.0 ); // initialise pair anisotropy constants to zero

               // set default uniaxial and cubic directions
               ku_vector.resize( 3 ); // resize to three elements

               ku_vector[ 0 ] = 0.0; // set direction along [0,0,1]
               ku_vector[ 1 ] = 0.0;
               ku_vector[ 2 ] = 1.0;

               // set default axis from which rotational angle phi is measured
               kr_vector.resize( 3 ); // resize to three elements
               
               kr_vector[ 0 ] = 1.0; // set direction along [1,0,0]
               kr_vector[ 1 ] = 0.0;
               kr_vector[ 2 ] = 0.0;

               // set default (y) axis perpendicular to both ku_vector and kr_vector
               kl_vector.resize( 3 );

               kl_vector[ 0 ] = 0.0;
               kl_vector[ 1 ] = 1.0;
               kl_vector[ 2 ] = 0.0;

               const double oneosqrt2 = 1.0 / sqrt( 2.0 );

               u1_vector.resize( 3 );  // resize to three elements
               u2_vector.resize( 3 );  // resize to three elements

               u1_vector[ 0 ] = oneosqrt2 * 1.0; // set direction along [1,1,0]
               u1_vector[ 1 ] = oneosqrt2 * 1.0;
               u1_vector[ 2 ] = 0.0;      

               u2_vector[ 0 ] = oneosqrt2 * 1.0;   // set direction along [1,-1,0]
               u2_vector[ 1 ] = - oneosqrt2 * 1.0;
               u2_vector[ 2 ] = 0.0;                             

               // set default uniaxial and cubic directions
               kc_vector1.resize( 3 ); // resize to three elements
               kc_vector2.resize( 3 ); // resize to three elements
               kc_vector3.resize( 3 ); // resize to three elements

               kc_vector1[ 0 ] = 1.0; // set direction alon [1,0,0]
               kc_vector1[ 1 ] = 0.0;
               kc_vector1[ 2 ] = 0.0;

               kc_vector2[ 0 ] = 0.0; // set direction alon [0,1,0]
               kc_vector2[ 1 ] = 1.0;
               kc_vector2[ 2 ] = 0.0;

               kc_vector3[ 0 ] = 0.0; // set direction alon [0,0,1]
               kc_vector3[ 1 ] = 0.0;
               kc_vector3[ 2 ] = 1.0;

               // set tensors as empty by default
               ku_tensor.resize( 9, 0.0 );
               kc_tensor.resize( 9, 0.0 );

            }; // end of constructor

      }; // end of anisotropy::internal::mp class

      //-----------------------------------------------------------------------------
      // Internal shared variables used for creation
      //-----------------------------------------------------------------------------

      extern std::vector<internal::mp_t> mp; // array of material properties

      extern bool initialised; // check module has been initialised

      // Flags for tesseral harmonic (or real spherical harmonic) expansion anisotorpy terms 
      extern bool enable_uniaxial_second_order;    // Flag to enable calculation of second order uniaxial anisotropy
      extern bool enable_rotational_2_1_order;     // Flag to enable calculation of second order theta first order phi anisotropy
      extern bool enable_rotational_2_1_order_odd; // Flag to enable calculation of second order theta first order phi odd anisotropy
      extern bool enable_rotational_2_2_order;     // Flag to enable calculation of second order theta second order phi anisotropy
      extern bool enable_rotational_2_2_order_odd; // Flag to enable calculation of second order theta second order phi odd anisotropy
      extern bool enable_uniaxial_fourth_order;    // Flag to enable calculation of fourth order uniaxial anisotropy
      extern bool enable_rotational_4_1_order;     // Flag to enable calculation of fourth order theta first order phi anisotropy
      extern bool enable_rotational_4_1_order_odd; // Flag to enable calculation of fourth order theta first order phi odd anisotropy
      extern bool enable_rotational_4_2_order;     // Flag to enable calculation of fourth order theta second order phi anisotropy
      extern bool enable_rotational_4_2_order_odd; // Flag to enable calculation of fourth order theta second order phi odd anisotropy
      extern bool enable_rotational_4_3_order;     // Flag to enable calculation of fourth order theta third order phi anisotropy
      extern bool enable_rotational_4_3_order_odd; // Flag to enable calculation of fourth order theta third order phi odd anisotropy
      extern bool enable_rotational_4_4_order;     // Flag to enable calculation of fourth order theta fourth order phi anisotropy
      extern bool enable_rotational_4_4_order_odd; // Flag to enable calculation of fourth order theta fourth order phi odd anisotropy
      extern bool enable_uniaxial_sixth_order;     // Flag to enable calculation of sixth order uniaxial anisotropy
      extern bool enable_rotational_6_1_order;     // Flag to enable calculation of sixth order theta first order phi anisotropy
      extern bool enable_rotational_6_1_order_odd; // Flag to enable caluclation of sixth order theta first order phi odd anisotropy
      extern bool enable_rotational_6_2_order;     // Flag to enable calculation of sixth order theta second order phi anisotropy
      extern bool enable_rotational_6_2_order_odd; // Flag to enable calculation of sixth order theta second order phi odd anisotropy
      extern bool enable_rotational_6_3_order;     // Flag to enable caluclation of sixth order theta third order phi anisotropy
      extern bool enable_rotational_6_3_order_odd; // Flag to enable calculation of sixth order theta third order phi odd anisotropy
      extern bool enable_rotational_6_4_order;     // Flag to enable calculation of sixth order theta fourth order phi anisotropy
      extern bool enable_rotational_6_4_order_odd; // Flag to enable calculation of sixth order theta fourth order phi odd anisotropy
      extern bool enable_rotational_6_5_order;     // Flag to enable calculation of sixth order theta fifth order phi anisotropy
      extern bool enable_rotational_6_5_order_odd; // Flag to enable calculation of sixth order theta fifth order phi odd anisotropy
      extern bool enable_rotational_6_6_order;     // Flag to enable calculation of sixth order theta fourth order phi anisotropy
      extern bool enable_rotational_6_6_order_odd; // Flag to enable calculation of sixth order theta fourth order phi odd anisotropy

      // Flags for biaxial anisotropy terms
      extern bool enable_biaxial_fourth_order_simple; // Flag to enable calculation of simplified fourth order biaxial anisotropy

      // Flags for cubic anisotropy terms
      extern bool enable_cubic_fourth_order;    // Flag to enable calculation of fourth order cubic anisotropy
      extern bool enable_cubic_sixth_order;     // Flag to enable calculation of sixth order cubic  anisotropy
      extern bool enable_cubic_fourth_order_rotation; // Flag to enable calculation of rotated cubic anisotropy

      // Flags for triaxial anisotropy
      extern bool enable_triaxial_anisotropy_rotated;
      extern bool enable_triaxial_fourth_order_rotated;
      extern bool enable_triaxial_anisotropy;
      extern bool enable_triaxial_fourth_order;

      // Flags for other anisotropies
      extern bool enable_neel_anisotropy; // Flag to turn on Neel anisotropy calculation (memory intensive at startup)
      extern bool enable_lattice_anisotropy; // Flag to turn on lattice anisotropy calculation
      extern bool enable_random_anisotropy; // Flag to enable random anisitropy initialisation

      // arrays for storing 1D collapsed Neel tensor
      extern std::vector< double > neel_tensor;

      // arrays for storing unrolled spherical-harmonic anisotropy constants in Tesla
      extern std::vector< double > ku2;
      extern std::vector< double > k2r1;
      extern std::vector< double > k2r1_odd;
      extern std::vector< double > k2r2;
      extern std::vector< double > k2r2_odd;
      extern std::vector< double > ku4;
      extern std::vector< double > k4r1;
      extern std::vector< double > k4r1_odd;
      extern std::vector< double > k4r2;
      extern std::vector< double > k4r2_odd;
      extern std::vector< double > k4r3;
      extern std::vector< double > k4r3_odd;
      extern std::vector< double > k4r4;
      extern std::vector< double > k4r4_odd;
      extern std::vector< double > ku6;
      extern std::vector< double > k6r1;
      extern std::vector< double > k6r1_odd;
      extern std::vector< double > k6r2;
      extern std::vector< double > k6r2_odd;
      extern std::vector< double > k6r3;
      extern std::vector< double > k6r3_odd;
      extern std::vector< double > k6r4;
      extern std::vector< double > k6r4_odd;
      extern std::vector< double > k6r5;
      extern std::vector< double > k6r5_odd;
      extern std::vector< double > k6r6;
      extern std::vector< double > k6r6_odd;

      extern std::vector< double > kc4;
      extern std::vector< double > kc6;

      // unrolled arrays for storing easy axes for each material
      extern std::vector< evec_t > ku_vector; // 001 easy axis direction
      extern std::vector< evec_t > kr_vector; // 100 phi = 0 direction
      extern std::vector< evec_t > kl_vector; // 010 phi = 90 direction

      extern bool native_neel_anisotropy_threshold;  // enables site-dependent surface threshold
      extern unsigned int neel_anisotropy_threshold; // global threshold for surface atoms
      extern double nearest_neighbour_distance;      // Control surface anisotropy nearest neighbour distance
      extern bool neel_range_dependent;              // Enable range dependent Neel anisotropy Lij = L0 exp(-F(r-r0)/r0)
      extern double neel_exponential_range;          // r0 value for range dependence of Neel anisotropy
      extern double neel_exponential_factor;         // F value for range dependence of Neel anisotropy

      extern std::vector< bool > triaxial_second_order_fixed_basis;
      extern std::vector< bool > triaxial_fourth_order_fixed_basis;

      extern std::vector< double > ku_triaxial_vector_x; // unit vector defining axis for uniaxial anisotropy
      extern std::vector< double > ku_triaxial_vector_y; // unit vector defining axis for uniaxial anisotropy
      extern std::vector< double > ku_triaxial_vector_z; // unit vector defining axis for uniaxial anisotropy

      extern std::vector< double > ku4_triaxial_vector_x; // unit vector defining axis for uniaxial anisotropy
      extern std::vector< double > ku4_triaxial_vector_y; // unit vector defining axis for uniaxial anisotropy
      extern std::vector< double > ku4_triaxial_vector_z; // unit vector defining axis for uniaxial anisotropy

      //basis vectors for second order triaxial - must be orthogonality
      extern std::vector < double > ku_triaxial_basis1x;
      extern std::vector < double > ku_triaxial_basis1y;
      extern std::vector < double > ku_triaxial_basis1z;

      extern std::vector < double > ku_triaxial_basis2x;
      extern std::vector < double > ku_triaxial_basis2y;
      extern std::vector < double > ku_triaxial_basis2z;

      extern std::vector < double > ku_triaxial_basis3x;
      extern std::vector < double > ku_triaxial_basis3y;
      extern std::vector < double > ku_triaxial_basis3z;

      //basis vectors for fourth order triaxial - must be orthogonality
      extern std::vector < double > ku4_triaxial_basis1x;
      extern std::vector < double > ku4_triaxial_basis1y;
      extern std::vector < double > ku4_triaxial_basis1z;

      extern std::vector < double > ku4_triaxial_basis2x;
      extern std::vector < double > ku4_triaxial_basis2y;
      extern std::vector < double > ku4_triaxial_basis2z;

      extern std::vector < double > ku4_triaxial_basis3x;
      extern std::vector < double > ku4_triaxial_basis3y;
      extern std::vector < double > ku4_triaxial_basis3z;

      // arrays for storing unrolled parameters for lattice anisotropy
      extern std::vector< double > klattice_array; // anisoptropy constant

      //-------------------------------------------------------------------------
      // internal function declarations
      //-------------------------------------------------------------------------
      
      // Fields
      void uniaxial_second_order_fields(  std::vector<double>& spin_array_x,
                                          std::vector<double>& spin_array_y,
                                          std::vector<double>& spin_array_z,
                                          std::vector<int>&    atom_material_array,
                                          std::vector<double>& field_array_x,
                                          std::vector<double>& field_array_y,
                                          std::vector<double>& field_array_z,
                                          const int start_index,
                                          const int end_index );

      void second_order_theta_first_order_phi_fields( std::vector<double>& spin_array_x,
                                                      std::vector<double>& spin_array_y,
                                                      std::vector<double>& spin_array_z,
                                                      std::vector<int>&    atom_material_array,
                                                      std::vector<double>& field_array_x,
                                                      std::vector<double>& field_array_y,
                                                      std::vector<double>& field_array_z,
                                                      const int start_index,
                                                      const int end_index );


      void second_order_theta_first_order_phi_odd_fields(   std::vector<double>& spin_array_x,
                                                            std::vector<double>& spin_array_y,
                                                            std::vector<double>& spin_array_z,
                                                            std::vector<int>&    atom_material_array,
                                                            std::vector<double>& field_array_x,
                                                            std::vector<double>& field_array_y,
                                                            std::vector<double>& field_array_z,
                                                            const int start_index,
                                                            const int end_index );

      void second_order_theta_second_order_phi_fields(   std::vector<double>& spin_array_x,
                                                         std::vector<double>& spin_array_y,
                                                         std::vector<double>& spin_array_z,
                                                         std::vector<int>&    atom_material_array,
                                                         std::vector<double>& field_array_x,
                                                         std::vector<double>& field_array_y,
                                                         std::vector<double>& field_array_z,
                                                         const int start_index,
                                                         const int end_index );

      void second_order_theta_second_order_phi_odd_fields(std::vector<double>& spin_array_x,
                                                      std::vector<double>& spin_array_y,
                                                      std::vector<double>& spin_array_z,
                                                      std::vector<int>&    atom_material_array,
                                                      std::vector<double>& field_array_x,
                                                      std::vector<double>& field_array_y,
                                                      std::vector<double>& field_array_z,
                                                      const int start_index,
                                                      const int end_index );

      void uniaxial_fourth_order_fields(  std::vector<double>& spin_array_x,
                                          std::vector<double>& spin_array_y,
                                          std::vector<double>& spin_array_z,
                                          std::vector<int>&    atom_material_array,
                                          std::vector<double>& field_array_x,
                                          std::vector<double>& field_array_y,
                                          std::vector<double>& field_array_z,
                                          const int start_index,
                                          const int end_index );

      void fourth_order_theta_first_order_phi_fields( std::vector<double>& spin_array_x,
                                                      std::vector<double>& spin_array_y,
                                                      std::vector<double>& spin_array_z,
                                                      std::vector<int>&    atom_material_array,
                                                      std::vector<double>& field_array_x,
                                                      std::vector<double>& field_array_y,
                                                      std::vector<double>& field_array_z,
                                                      const int start_index,
                                                      const int end_index );

      void fourth_order_theta_first_order_phi_odd_fields(   std::vector<double>& spin_array_x,
                                                            std::vector<double>& spin_array_y,
                                                            std::vector<double>& spin_array_z,
                                                            std::vector<int>&    atom_material_array,
                                                            std::vector<double>& field_array_x,
                                                            std::vector<double>& field_array_y,
                                                            std::vector<double>& field_array_z,
                                                            const int start_index,
                                                            const int end_index );

      void fourth_order_theta_second_order_phi_fields(   std::vector<double>& spin_array_x,
                                                         std::vector<double>& spin_array_y,
                                                         std::vector<double>& spin_array_z,
                                                         std::vector<int>&    atom_material_array,
                                                         std::vector<double>& field_array_x,
                                                         std::vector<double>& field_array_y,
                                                         std::vector<double>& field_array_z,
                                                         const int start_index,
                                                         const int end_index );

      void fourth_order_theta_second_order_phi_odd_fields(  std::vector<double>& spin_array_x,
                                                            std::vector<double>& spin_array_y,
                                                            std::vector<double>& spin_array_z,
                                                            std::vector<int>&    atom_material_array,
                                                            std::vector<double>& field_array_x,
                                                            std::vector<double>& field_array_y,
                                                            std::vector<double>& field_array_z,
                                                            const int start_index,
                                                            const int end_index );

      void fourth_order_theta_third_order_phi_fields( std::vector<double>& spin_array_x,
                                                      std::vector<double>& spin_array_y,
                                                      std::vector<double>& spin_array_z,
                                                      std::vector<int>&    atom_material_array,
                                                      std::vector<double>& field_array_x,
                                                      std::vector<double>& field_array_y,
                                                      std::vector<double>& field_array_z,
                                                      const int start_index,
                                                      const int end_index );

      void fourth_order_theta_third_order_phi_odd_fields(   std::vector<double>& spin_array_x,
                                                            std::vector<double>& spin_array_y,
                                                            std::vector<double>& spin_array_z,
                                                            std::vector<int>&    atom_material_array,
                                                            std::vector<double>& field_array_x,
                                                            std::vector<double>& field_array_y,
                                                            std::vector<double>& field_array_z,
                                                            const int start_index,
                                                            const int end_index );

      void fourth_order_theta_fourth_order_phi_fields(   std::vector<double>& spin_array_x,
                                                         std::vector<double>& spin_array_y,
                                                         std::vector<double>& spin_array_z,
                                                         std::vector<int>&    atom_material_array,
                                                         std::vector<double>& field_array_x,
                                                         std::vector<double>& field_array_y,
                                                         std::vector<double>& field_array_z,
                                                         const int start_index,
                                                         const int end_index );

      void fourth_order_theta_fourth_order_phi_odd_fields(  std::vector<double>& spin_array_x,
                                                            std::vector<double>& spin_array_y,
                                                            std::vector<double>& spin_array_z,
                                                            std::vector<int>&    atom_material_array,
                                                            std::vector<double>& field_array_x,
                                                            std::vector<double>& field_array_y,
                                                            std::vector<double>& field_array_z,
                                                            const int start_index,
                                                            const int end_index );

      void uniaxial_sixth_order_fields(   std::vector<double>& spin_array_x,
                                          std::vector<double>& spin_array_y,
                                          std::vector<double>& spin_array_z,
                                          std::vector<int>&    atom_material_array,
                                          std::vector<double>& field_array_x,
                                          std::vector<double>& field_array_y,
                                          std::vector<double>& field_array_z,
                                          const int start_index,
                                          const int end_index );

      void sixth_order_theta_first_order_phi_fields(  std::vector<double>& spin_array_x,
                                                      std::vector<double>& spin_array_y,
                                                      std::vector<double>& spin_array_z,
                                                      std::vector<int>&    atom_material_array,
                                                      std::vector<double>& field_array_x,
                                                      std::vector<double>& field_array_y,
                                                      std::vector<double>& field_array_z,
                                                      const int start_index,
                                                      const int end_index );

      void sixth_order_theta_first_order_phi_odd_fields( std::vector<double>& spin_array_x,
                                                         std::vector<double>& spin_array_y,
                                                         std::vector<double>& spin_array_z,
                                                         std::vector<int>&    atom_material_array,
                                                         std::vector<double>& field_array_x,
                                                         std::vector<double>& field_array_y,
                                                         std::vector<double>& field_array_z,
                                                         const int start_index,
                                                         const int end_index );

      void sixth_order_theta_second_order_phi_fields( std::vector<double>& spin_array_x,
                                                      std::vector<double>& spin_array_y,
                                                      std::vector<double>& spin_array_z,
                                                      std::vector<int>&    atom_material_array,
                                                      std::vector<double>& field_array_x,
                                                      std::vector<double>& field_array_y,
                                                      std::vector<double>& field_array_z,
                                                      const int start_index,
                                                      const int end_index );

      void sixth_order_theta_second_order_phi_odd_fields(   std::vector<double>& spin_array_x,
                                                            std::vector<double>& spin_array_y,
                                                            std::vector<double>& spin_array_z,
                                                            std::vector<int>&    atom_material_array,
                                                            std::vector<double>& field_array_x,
                                                            std::vector<double>& field_array_y,
                                                            std::vector<double>& field_array_z,
                                                            const int start_index,
                                                            const int end_index );

      void sixth_order_theta_third_order_phi_fields(  std::vector<double>& spin_array_x,
                                                      std::vector<double>& spin_array_y,
                                                      std::vector<double>& spin_array_z,
                                                      std::vector<int>&    atom_material_array,
                                                      std::vector<double>& field_array_x,
                                                      std::vector<double>& field_array_y,
                                                      std::vector<double>& field_array_z,
                                                      const int start_index,
                                                      const int end_index );

      void sixth_order_theta_third_order_phi_odd_fields(   std::vector<double>& spin_array_x,
                                                            std::vector<double>& spin_array_y,
                                                            std::vector<double>& spin_array_z,
                                                            std::vector<int>&    atom_material_array,
                                                            std::vector<double>& field_array_x,
                                                            std::vector<double>& field_array_y,
                                                            std::vector<double>& field_array_z,
                                                            const int start_index,
                                                            const int end_index );

      void sixth_order_theta_fourth_order_phi_fields( std::vector<double>& spin_array_x,
                                                      std::vector<double>& spin_array_y,
                                                      std::vector<double>& spin_array_z,
                                                      std::vector<int>&    atom_material_array,
                                                      std::vector<double>& field_array_x,
                                                      std::vector<double>& field_array_y,
                                                      std::vector<double>& field_array_z,
                                                      const int start_index,
                                                      const int end_index );

      void sixth_order_theta_fourth_order_phi_odd_fields(   std::vector<double>& spin_array_x,
                                                            std::vector<double>& spin_array_y,
                                                            std::vector<double>& spin_array_z,
                                                            std::vector<int>&    atom_material_array,
                                                            std::vector<double>& field_array_x,
                                                            std::vector<double>& field_array_y,
                                                            std::vector<double>& field_array_z,
                                                            const int start_index,
                                                            const int end_index );

      void sixth_order_theta_fifth_order_phi_fields(  std::vector<double>& spin_array_x,
                                                      std::vector<double>& spin_array_y,
                                                      std::vector<double>& spin_array_z,
                                                      std::vector<int>&    atom_material_array,
                                                      std::vector<double>& field_array_x,
                                                      std::vector<double>& field_array_y,
                                                      std::vector<double>& field_array_z,
                                                      const int start_index,
                                                      const int end_index );

      void sixth_order_theta_fifth_order_phi_odd_fields( std::vector<double>& spin_array_x,
                                                         std::vector<double>& spin_array_y,
                                                         std::vector<double>& spin_array_z,
                                                         std::vector<int>&    atom_material_array,
                                                         std::vector<double>& field_array_x,
                                                         std::vector<double>& field_array_y,
                                                         std::vector<double>& field_array_z,
                                                         const int start_index,
                                                         const int end_index );

      void sixth_order_theta_sixth_order_phi_fields(  std::vector<double>& spin_array_x,
                                                      std::vector<double>& spin_array_y,
                                                      std::vector<double>& spin_array_z,
                                                      std::vector<int>&    atom_material_array,
                                                      std::vector<double>& field_array_x,
                                                      std::vector<double>& field_array_y,
                                                      std::vector<double>& field_array_z,
                                                      const int start_index,
                                                      const int end_index );

      void sixth_order_theta_sixth_order_phi_odd_fields( std::vector<double>& spin_array_x,
                                                         std::vector<double>& spin_array_y,
                                                         std::vector<double>& spin_array_z,
                                                         std::vector<int>&    atom_material_array,
                                                         std::vector<double>& field_array_x,
                                                         std::vector<double>& field_array_y,
                                                         std::vector<double>& field_array_z,
                                                         const int start_index,
                                                         const int end_index );

      void triaxial_second_order_fields_fixed_basis(std::vector<double>& spin_array_x,
                                                     std::vector<double>& spin_array_y,
                                                     std::vector<double>& spin_array_z,
                                                     std::vector<int>&    atom_material_array,
                                                     std::vector<double>& field_array_x,
                                                     std::vector<double>& field_array_y,
                                                     std::vector<double>& field_array_z,
                                                     const int start_index,
                                                     const int end_index);

      void triaxial_second_order_fields(std::vector<double>& spin_array_x,
                                                    std::vector<double>& spin_array_y,
                                                    std::vector<double>& spin_array_z,
                                                    std::vector<int>&    atom_material_array,
                                                    std::vector<double>& field_array_x,
                                                    std::vector<double>& field_array_y,
                                                    std::vector<double>& field_array_z,
                                                    const int start_index,
                                                    const int end_index);

      void triaxial_fourth_order_fields_fixed_basis(std::vector<double>& spin_array_x,
                                                  std::vector<double>& spin_array_y,
                                                  std::vector<double>& spin_array_z,
                                                  std::vector<int>&    atom_material_array,
                                                  std::vector<double>& field_array_x,
                                                  std::vector<double>& field_array_y,
                                                  std::vector<double>& field_array_z,
                                                  const int start_index,
                                                  const int end_index);

      void triaxial_fourth_order_fields(std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index);

      void biaxial_fourth_order_simple_fields(std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index);                                        

      void cubic_fourth_order_fields(std::vector<double>& spin_array_x,
                                     std::vector<double>& spin_array_y,
                                     std::vector<double>& spin_array_z,
                                     std::vector<int>&    atom_material_array,
                                     std::vector<double>& field_array_x,
                                     std::vector<double>& field_array_y,
                                     std::vector<double>& field_array_z,
                                     const int start_index,
                                     const int end_index);

      void cubic_fourth_order_rotation_fields(std::vector<double>& spin_array_x,
                                              std::vector<double>& spin_array_y,
                                              std::vector<double>& spin_array_z,
                                              std::vector<int>&    atom_material_array,
                                              std::vector<double>& field_array_x,
                                              std::vector<double>& field_array_y,
                                              std::vector<double>& field_array_z,
                                              const int start_index,
                                              const int end_index);

      void cubic_sixth_order_fields( std::vector<double>& spin_array_x,
                                     std::vector<double>& spin_array_y,
                                     std::vector<double>& spin_array_z,
                                     std::vector<int>&    atom_material_array,
                                     std::vector<double>& field_array_x,
                                     std::vector<double>& field_array_y,
                                     std::vector<double>& field_array_z,
                                     const int start_index,
                                     const int end_index);

      void neel_fields( std::vector<double>& spin_array_x,
                        std::vector<double>& spin_array_y,
                        std::vector<double>& spin_array_z,
                        std::vector<int>&    atom_material_array,
                        std::vector<double>& field_array_x,
                        std::vector<double>& field_array_y,
                        std::vector<double>& field_array_z,
                        const int start_index,
                        const int end_index);

      void lattice_fields(std::vector<double>& spin_array_x,
                          std::vector<double>& spin_array_y,
                          std::vector<double>& spin_array_z,
                          std::vector<int>&    type_array,
                          std::vector<double>& field_array_x,
                          std::vector<double>& field_array_y,
                          std::vector<double>& field_array_z,
                          const int start_index,
                          const int end_index,
                          const double temperature);

      // Energies

      double uniaxial_second_order_energy(   const int atom,
                                             const int mat,
                                             const double sx,
                                             const double sy,
                                             const double sz );

      double second_order_theta_first_order_phi_energy(  const int atom,
                                                         const int mat,
                                                         const double sx,
                                                         const double sy,
                                                         const double sz );

      double second_order_theta_first_order_phi_odd_energy( const int atom,
                                                            const int mat,
                                                            const double sx,
                                                            const double sy,
                                                            const double sz );

      double second_order_theta_second_order_phi_energy( const int atom,
                                                         const int mat,
                                                         const double sx,
                                                         const double sy,
                                                         const double sz );

      double second_order_theta_second_order_phi_odd_energy(   const int atom,
                                                               const int mat,
                                                               const double sx,
                                                               const double sy,
                                                               const double sz );

      double uniaxial_fourth_order_energy(   const int atom,
                                             const int mat,
                                             const double sx,
                                             const double sy,
                                             const double sz );

      double fourth_order_theta_first_order_phi_energy(  const int atom,
                                                         const int mat,
                                                         const double sx,
                                                         const double sy,
                                                         const double sz );

      double fourth_order_theta_first_order_phi_odd_energy( const int atom,
                                                            const int mat,
                                                            const double sx,
                                                            const double sy,
                                                            const double sz );

      double fourth_order_theta_second_order_phi_energy( const int atom,
                                                         const int mat,
                                                         const double sx,
                                                         const double sy,
                                                         const double sz );

      double fourth_order_theta_second_order_phi_odd_energy(   const int atom,
                                                               const int mat,
                                                               const double sx,
                                                               const double sy,
                                                               const double sz );

      double fourth_order_theta_third_order_phi_odd_energy( const int atom,
                                                            const int mat,
                                                            const double sx,
                                                            const double sy,
                                                            const double sz );

      double fourth_order_theta_third_order_phi_energy(  const int atom,
                                                         const int mat,
                                                         const double sx,
                                                         const double sy,
                                                         const double sz );

      double fourth_order_theta_fourth_order_phi_energy( const int atom,
                                                         const int mat,
                                                         const double sx,
                                                         const double sy,
                                                         const double sz );

      double fourth_order_theta_fourth_order_phi_odd_energy(   const int atom,
                                                               const int mat,
                                                               const double sx,
                                                               const double sy,
                                                               const double sz );

      double uniaxial_sixth_order_energy( const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz);

      double sixth_order_theta_first_order_phi_energy(   const int atom,
                                                         const int mat,
                                                         const double sx,
                                                         const double sy,
                                                         const double sz );

      double sixth_order_theta_first_order_phi_odd_energy(  const int atom,
                                                            const int mat,
                                                            const double sx,
                                                            const double sy,
                                                            const double sz );

      double sixth_order_theta_second_order_phi_energy(  const int atom,
                                                         const int mat,
                                                         const double sx,
                                                         const double sy,
                                                         const double sz );

      double sixth_order_theta_second_order_phi_odd_energy( const int atom,
                                                            const int mat,
                                                            const double sx,
                                                            const double sy,
                                                            const double sz );

      double sixth_order_theta_third_order_phi_energy(   const int atom,
                                                         const int mat,
                                                         const double sx,
                                                         const double sy,
                                                         const double sz );

      double sixth_order_theta_third_order_phi_odd_energy(  const int atom,
                                                            const int mat,
                                                            const double sx,
                                                            const double sy,
                                                            const double sz );

      double sixth_order_theta_fourth_order_phi_energy(  const int atom,
                                                         const int mat,
                                                         const double sx,
                                                         const double sy,
                                                         const double sz );

      double sixth_order_theta_fourth_order_phi_odd_energy( const int atom,
                                                            const int mat,
                                                            const double sx,
                                                            const double sy,
                                                            const double sz );

      double sixth_order_theta_fifth_order_phi_energy(   const int atom,
                                                         const int mat,
                                                         const double sx,
                                                         const double sy,
                                                         const double sz );

      double sixth_order_theta_fifth_order_phi_odd_energy(  const int atom,
                                                            const int mat,
                                                            const double sx,
                                                            const double sy,
                                                            const double sz );

      double sixth_order_theta_sixth_order_phi_energy(   const int atom,
                                                         const int mat,
                                                         const double sx,
                                                         const double sy,
                                                         const double sz );

      double sixth_order_theta_sixth_order_phi_odd_energy(  const int atom,
                                                            const int mat,
                                                            const double sx,
                                                            const double sy,
                                                            const double sz );

      double triaxial_second_order_energy_fixed_basis(const int atom,
                                                      const int mat,
                                                      const double sx,
                                                      const double sy,
                                                      const double sz);

      double triaxial_fourth_order_energy_fixed_basis(const int atom,
                                                      const int mat,
                                                      const double sx,
                                                      const double sy,
                                                      const double sz);

      double triaxial_second_order_energy(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz);

      double triaxial_fourth_order_energy(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz);

      double biaxial_fourth_order_simple_energy(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz);                                          

      double cubic_fourth_order_energy(const int atom,
                                       const int mat,
                                       const double sx,
                                       const double sy,
                                       const double sz);

      double cubic_fourth_order_rotation_energy(const int atom,
                                       const int mat,
                                       const double sx,
                                       const double sy,
                                       const double sz);


      double cubic_sixth_order_energy( const int atom,
                                       const int mat,
                                       const double sx,
                                       const double sy,
                                       const double sz);

      double neel_energy( const int atom,
                          const int mat,
                          const double sx,
                          const double sy,
                          const double sz);

      double lattice_energy(const int atom, const int mat, const double sx, const double sy, const double sz, const double temperature);

      void initialise_neel_anisotropy_tensor(std::vector <std::vector <bool> >& nearest_neighbour_interactions_list,
                                             std::vector<std::vector <neighbours::neighbour_t> >& cneighbourlist);

   } // end of internal namespace

   //-------------------------------------------------------------------------
   // function declarations
   //-------------------------------------------------------------------------

   //-------------------------------------------------------------------------
   // simple inline function to convert atom,i,j into 1D tensor coordinates
   //-------------------------------------------------------------------------
   inline unsigned int index( const unsigned int atom, const unsigned int i, const unsigned int j )
   {

      return 9 * atom + 3 * i + j;

   }

} // end of anisotropy namespace

#endif //ANISOTROPY_INTERNAL_H_
