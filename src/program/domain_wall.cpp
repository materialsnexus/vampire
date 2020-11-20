//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and Andrea Meo 2014-2018. All rights reserved.
//
//-----------------------------------------------------------------------------
//

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

namespace program{

//------------------------------------------------------------------------------
// Program to calculate a simple time series
//------------------------------------------------------------------------------
void domain_wall(){

	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "program::domain walls has been called" << std::endl;

	double temp=sim::temperature;
	#ifdef MPICF
		int num_local_atoms =  vmpi::num_core_atoms+vmpi::num_bdry_atoms;
	#else
		int num_local_atoms = atoms::num_atoms;
	#endif

    int num_averages = 50;
    num_averages = num_averages/2.0;
    int num_dw_cells = (cs::system_dimensions[sim::domain_wall_axis])/sim::domain_wall_discretisation + 1;
    std::vector <double > atom_to_cell_array(num_local_atoms,0);
    std::vector < std::vector < double > > mag;
    std::vector < std::vector < double > > initial_mag;
    std::vector < std::vector < double > > counter;
    std::vector < std::vector < double > > sum_mag;
    std::vector < std::vector < double > > num_atoms_in_cell;

    sum_mag.resize(num_dw_cells*3);
    num_atoms_in_cell.resize(num_dw_cells*3);
    counter.resize(num_dw_cells*3);
    initial_mag.resize(num_dw_cells*3);
    mag.resize(num_dw_cells*3);

    for (int i = 0; i < num_dw_cells*3; ++i)
    {
      sum_mag[i].resize(mp::num_materials);
       num_atoms_in_cell[i].resize(mp::num_materials);
       counter[i].resize(mp::num_materials);
       initial_mag[i].resize(mp::num_materials);
       mag[i].resize(mp::num_materials);
    }

   //works out which atoms are in which cells and sets cells based on wether the domain wall
   //is along x or y or z

   if (sim::domain_wall_axis == 0){
      for(int atom=0;atom<num_local_atoms;atom++){
         int mat = atoms::type_array[atom];
         int cell = atoms::x_coord_array[atom]/sim::domain_wall_discretisation;
         atom_to_cell_array[atom] = cell;
         num_atoms_in_cell[cell][mat] ++;
		//	if (cell > num_dw_cells) std::cout << atoms::x_coord_array[atom] << '\t' << sim::domain_wall_discretisation <<std::endl;
      }
   }
   if (sim::domain_wall_axis == 1){
      for(int atom=0;atom<num_local_atoms;atom++){
         int cell = atoms::y_coord_array[atom]/sim::domain_wall_discretisation;
         atom_to_cell_array[atom] = cell;

      }
   }
   if (sim::domain_wall_axis == 2){
      for(int atom=0;atom<num_local_atoms;atom++){
         int cell = atoms::z_coord_array[atom]/sim::domain_wall_discretisation;
         atom_to_cell_array[atom] = cell;
      }
   }


   //reverses the magentisation of atoms further away than the domain wall distance.
   if (!sim::load_checkpoint_flag){
      if (sim::domain_wall_axis == 0){
         for(int atom=0;atom<num_local_atoms;atom++){
            if (atoms::x_coord_array[atom] > cs::system_dimensions[0]*sim::domain_wall_position){
               atoms::x_spin_array[atom] = -atoms::x_spin_array[atom];
               atoms::y_spin_array[atom] = -atoms::y_spin_array[atom];
      	      atoms::z_spin_array[atom] = -atoms::z_spin_array[atom];
             }
         }
      }

      if (sim::domain_wall_axis == 1){
         for(int atom=0;atom<num_local_atoms;atom++){
            if (atoms::y_coord_array[atom] > cs::system_dimensions[1]*sim::domain_wall_position){
               atoms::x_spin_array[atom] = -atoms::x_spin_array[atom];
               atoms::y_spin_array[atom] = -atoms::y_spin_array[atom];
      	      atoms::z_spin_array[atom] = -atoms::z_spin_array[atom];
             }
         }
      }

      if (sim::domain_wall_axis == 2){
         for(int atom=0;atom<num_local_atoms;atom++){
            if (atoms::z_coord_array[atom] > cs::system_dimensions[2]*sim::domain_wall_position){
               atoms::x_spin_array[atom] = -atoms::x_spin_array[atom];
               atoms::y_spin_array[atom] = -atoms::y_spin_array[atom];
      	      atoms::z_spin_array[atom] = -atoms::z_spin_array[atom];
             }
         }
      }
   }

	if (sim::anti_PBC[0] || sim::anti_PBC[1] || sim::anti_PBC[2]){
		for (int atom = 0; atom <num_local_atoms; atom++){
			const int start = atoms::neighbour_list_start_index[atom];
			const int end   = atoms::neighbour_list_end_index[atom] + 1;
			for(int nn=start;nn<end;nn++){
				const int natom = atoms::neighbour_list_array[nn];
				if (sim::anti_PBC[0] == true){
					const double dx = atoms::x_coord_array[atom] - atoms::x_coord_array[natom];
					if (dx > cs::system_dimensions[0]-5){
						atoms::i_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij = -1.0*atoms::i_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij;
					}
				}
				if (sim::anti_PBC[1] == true){

					const double dy = atoms::y_coord_array[atom] - atoms::y_coord_array[natom];
					if (dy > cs::system_dimensions[1]-20){
					atoms::i_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij = -1.0*atoms::i_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij;
					}
				}
				if (sim::anti_PBC[2] == true){
					const double dz = atoms::z_coord_array[atom] - atoms::z_coord_array[natom];
					if (dz > cs::system_dimensions[2]-20){
					atoms::i_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij = -1.0*atoms::i_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij * -1;
				}
				}
			}
		}
	}

//std::cout << "end loop" << std::endl;

   // Set equilibration temperature only if continue checkpoint not loaded
   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
   else{
	   // Set equilibration temperature
	   sim::temperature=sim::Teq;
   }

	// Equilibrate system
	while(sim::time<sim::equilibration_time){

		sim::integrate(sim::partial_time);

		// Calculate magnetisation statistics
		stats::mag_m();

		// Output data
		vout::data();
	}

///	std::cout << "end equilbration" << std::endl;


   // Set temperature and reset stats only if continue checkpoint not loaded
   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
   else{

      // set simulation temperature
	   sim::temperature = temp;

      // Reset mean magnetisation counters
      stats::mag_m_reset();

   }

//	std::cout << "reset mag" << std::endl;


	// Perform Time Series
	while(sim::time<sim::equilibration_time+sim::total_time){

      for (int cell = 0; cell < num_dw_cells; cell++){
         for (int mat = 0; mat < mp::num_materials; mat ++){
      //      std::cout << mat << '\t' << cell << "\t" << mag[cell*3 + 0][mat] << "\t" << mag[cell*3 + 1][mat] << "\t" << mag[cell*3 + 2][mat] << std::endl;
            mag[3*cell + 0][mat] = 0;
            mag[3*cell + 1][mat] = 0;
            mag[3*cell + 2][mat] = 0;
            initial_mag[3*cell + 0][mat] = 0;
            initial_mag[3*cell + 1][mat] = 0;
            initial_mag[3*cell + 2][mat] = 0;
         }
      }

	//	std::cout << "initially set cells" << std::endl;

		// Integrate system
		sim::integrate(sim::partial_time);

      // Calculate magnetisation statistics
		stats::mag_m();
//		std::cout << "integrate" <<  "\t" << num_local_atoms << std::endl;


      std::ofstream myfile;
      string filename = "dw-" + std::to_string(sim::time) + ".txt";
      myfile.open (filename);

      for(int atom=0;atom<num_local_atoms;atom++){
         int cell = atom_to_cell_array[atom];
         int mat = atoms::type_array[atom];
		//	std::cout << atom << '\t' << cell << "\t" << mat << "\t" << num_local_atoms << "\t" << atoms::x_spin_array.size() << "\t" << initial_mag.size() << "\t" << initial_mag[0].size() << std::endl;
         initial_mag[cell*3 + 0][mat] += atoms::x_spin_array[atom];
	//		std::cout << "a" << std::endl;
         initial_mag[cell*3 + 1][mat] += atoms::y_spin_array[atom];
	//		std::cout << "a" << std::endl;

         initial_mag[cell*3 + 2][mat] += atoms::z_spin_array[atom];
	//		std::cout << "a" << std::endl;

      //   if (mat == 0) myfile << cell << "\t" << mag[cell*3 + 0][mat] << "\t" << mag[cell*3 + 1][mat] << "\t" << mag[cell*3 + 2][mat] << std::endl;

      }
	//	std::cout << "loop over atoms" << std::endl;

	//	for (int cell = 0; cell < num_dw_cells; cell++){
       for (int cell = num_averages; cell < num_dw_cells - num_averages; cell++){
          int mat = 0;
         // mag[cell*3 + 0][mat] = initial_mag[cell*3 + 0][mat]/num_atoms_in_cell[cell][mat];
         // mag[cell*3 + 1][mat] = initial_mag[cell*3 + 1][mat]/num_atoms_in_cell[cell][mat];
         // mag[cell*3 + 2][mat] = initial_mag[cell*3 + 2][mat]/num_atoms_in_cell[cell][mat];
         for (int i = cell - num_averages; i < cell + num_averages; i ++  ){
            mag[cell*3 + 0][mat] = mag[cell*3 + 0][mat] + initial_mag[i*3 + 0][mat]/(num_atoms_in_cell[i][mat]*num_averages*2+1);
            mag[cell*3 + 1][mat] = mag[cell*3 + 1][mat] + initial_mag[i*3 + 1][mat]/(num_atoms_in_cell[i][mat]*num_averages*2+1);
            mag[cell*3 + 2][mat] = mag[cell*3 + 2][mat] + initial_mag[i*3 + 2][mat]/(num_atoms_in_cell[i][mat]*num_averages*2+1);
          }
         //mag[cell*3 + 0][mat]= mag[cell*3 + 0][mat]/(num_averages+1);
         myfile << cell << '\t'<< '\t' << mag[cell*3 + 0][mat] << "\t" << mag[cell*3 + 1][mat] << "\t" << mag[cell*3 + 2][mat] <<"\t" << std::endl;// av_dl << std::endl;//'\t' << sum_mag[new_pos*3 + 0][mat] << "\t" << sum_mag[new_pos*3 + 1][mat] << "\t" << sum_mag[new_pos*3 + 2][mat] << "\t" << counter[new_pos][mat] << "\t" << dl[cell][mat] << std::endl;

      }
		// for (int cell = num_averages +1; cell < num_dw_cells - num_averages ; cell++){
		// 	if (sim::anti_PBC[0] == true){
		// 		val = mag[cell*3 + 0][0]*mag[(cell-1)*3 + 0][0];
		// 	}
		// }
		//
      // double max_dl = 0;
      // double centre = 0;
      // std::vector < std::vector < double > > dl;
      // dl.resize(num_dw_cells);
      // for (int i = 0; i < num_dw_cells; ++i)
      //     dl[i].resize(mp::num_materials);
		//
      // for (int cell = num_averages +1; cell < num_dw_cells - num_averages ; cell++){
      //       int mat = 0;
      //       double dx = mag[cell*3 + 0][0] -  mag[(cell-1)*3 + 0][0];
      //       double dy = mag[cell*3 + 1][0] -  mag[(cell-1)*3 + 1][0];
      //       double dz = mag[cell*3 + 2][0] -  mag[(cell-1)*3 + 2][0];
      //       dl[cell][mat] = dx*dx + dy*dy + dz*dz;
      // }
		//
      // for (int cell = 2*num_averages + 1; cell < num_dw_cells - 2*num_averages; cell++){
      //    int mat = 0;
      //    double av_dl = dl[cell][mat];
		//
      //    for (int i = cell - num_averages; i < cell + num_averages; i ++  ){
      //    //   std::cout << cell << '\t' << i << "\t" << dl[i][mat] << std::endl;
      //       av_dl = av_dl + dl[i][mat];
      //    }
      //    av_dl = av_dl/(num_averages*2 + 1);
      // //   std::cout << av_dl << std::endl;
      //     if (av_dl > max_dl) {
      //        max_dl = av_dl;
      //        centre = cell;
      //     }
      //    // myfile << cell << '\t'<< '\t' << mag[cell*3 + 0][mat] << "\t" << mag[cell*3 + 1][mat] << "\t" << mag[cell*3 + 2][mat] <<"\t" <<  av_dl << std::endl;//'\t' << sum_mag[new_pos*3 + 0][mat] << "\t" << sum_mag[new_pos*3 + 1][mat] << "\t" << sum_mag[new_pos*3 + 2][mat] << "\t" << counter[new_pos][mat] << "\t" << dl[cell][mat] << std::endl;
      //     }
      // std::cout << centre<< "\t" << sim::time << std::endl;
		//
		//
      // double real_centre = cs::system_dimensions[sim::domain_wall_axis]*0.5/sim::domain_wall_discretisation;
      // double dcentre = real_centre - centre;
      // //std::cout << dcentre << std::endl;
      // for (int cell = num_averages + 1; cell < num_dw_cells - num_averages; cell++){
      //    double new_pos = cell + dcentre;
      // //   std::cout << cell << '\t' << dcentre << '\'t << new_pos << std::endl;
      //    int mat = 0;
      // //   std::cout << cell << '\t' << centre << '\t' << real_centre << '\t' << dcentre << '\t' << new_pos << std::endl;
      //    if (new_pos > -0.01 && new_pos < num_dw_cells){
		//
      //       sum_mag[new_pos*3 + 0][mat] += mag[cell*3 + 0][mat];
      //       sum_mag[new_pos*3 + 1][mat] += mag[cell*3 + 1][mat];
      //       sum_mag[new_pos*3 + 2][mat] += mag[cell*3 + 2][mat];
      //       counter[new_pos][mat] ++;
      //       myfile << cell << '\t' << new_pos << '\t' << mag[cell*3 + 0][mat] << "\t" << mag[cell*3 + 1][mat] << "\t" << mag[cell*3 + 2][mat] << '\t' << sum_mag[new_pos*3 + 0][mat] << "\t" << sum_mag[new_pos*3 + 1][mat] << "\t" << sum_mag[new_pos*3 + 2][mat] << "\t" << counter[new_pos][mat] << "\t" << dl[cell][mat] << std::endl;
      //    }
      // }
      // sim::domain_wall_centre = centre*sim::domain_wall_discretisation/10;

		// Output data
		vout::data();
      myfile.close();
   //   std::cout << "a" <<std::endl;
	}

   //   std::cout << num_averages + 1 << "\t" << num_dw_cells - num_averages <<std::endl;
   // for (int cell = num_averages + 1; cell < num_dw_cells - num_averages; cell++){
   //    int mat = 0;
   //    if (counter[cell][mat] > 0){
   //       std::cout << cell*sim::domain_wall_discretisation << '\t' << sum_mag[cell*3 + 0][mat]/(counter[cell][mat]) << "\t" << sum_mag[cell*3 + 1][mat]/(counter[cell][mat]) << "\t" << sum_mag[cell*3 + 2][mat]/(counter[cell][mat])<< "\t" << counter[cell][mat] << std::endl;
   //    }
   // }

}

}//end of namespace program
