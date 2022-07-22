/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
#include <math.h>

void create_cell_types( void )
{
	// set the random seed
	SeedRandom( parameters.ints("random_seed") );

	/*
	   Put any modifications to default cell definition here if you
	   want to have "inherited" by other cell types.

	   This is a good place to set default functions.
	*/

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL;
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based;
	cell_defaults.functions.custom_cell_rule = NULL;
	cell_defaults.functions.contact_function = NULL;

	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;

	/*
	   This parses the cell definitions in the XML config file.
	*/

	initialize_cell_definitions_from_pugixml();

	/*
	   This builds the map of cell definitions and summarizes the setup.
	*/

	build_cell_definitions_maps();

	/*
	   This intializes cell signal and response dictionaries
	*/

	setup_signal_behavior_dictionaries();

	/*
	   Put any modifications to individual cell definitions here.

	   This is a good place to set custom functions.
	*/

	cell_defaults.functions.update_phenotype = phenotype_function;
	cell_defaults.functions.custom_cell_rule = custom_function;
	cell_defaults.functions.contact_function = contact_function;

	/*
	   This builds the map of cell definitions and summarizes the setup.
	*/

	display_cell_definitions( std::cout );

	return;
}

void setup_microenvironment( void )
{
	// set domain parameters

	// put any custom code to set non-homogeneous initial conditions or
	// extra Dirichlet nodes here.

	// initialize BioFVM

	initialize_microenvironment();

	return;
}

void setup_tissue_rwh( void )
{
	Cell* pC;
    Cell_Definition* pCD = cell_definitions_by_index[0];
    std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
        std::vector<double> position = {0,0,0};
        double xval = 1.0;
        double yval = 0.;
        double zval = 0.;
        position[0] = xval;
        position[1] = yval;
        position[2] = zval;
        std::cout << xval<<", "<< yval<<", "<< zval<<std::endl;
        pC = create_cell( *pCD );
        pC->assign_position( position );

    pCD = cell_definitions_by_index[1];
    std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
        xval = -1.0;
        yval = 0.;
        zval = 0.;
        position[0] = xval;
        position[1] = yval;
        position[2] = zval;
        std::cout << xval<<", "<< yval<<", "<< zval<<std::endl;
        pC = create_cell( *pCD );
        pC->assign_position( position );
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0];
	double Ymin = microenvironment.mesh.bounding_box[1];
	double Zmin = microenvironment.mesh.bounding_box[2];

	double Xmax = microenvironment.mesh.bounding_box[3];
	double Ymax = microenvironment.mesh.bounding_box[4];
	double Zmax = microenvironment.mesh.bounding_box[5];

	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0;
		Zmax = 0.0;
	}

	double Xrange = Xmax - Xmin;
	double Yrange = Ymax - Ymin;
	double Zrange = Zmax - Zmin;

	// create some of each type of cell

	Cell* pC;

    double cyl_end_xval = parameters.doubles("cyl_x1");  // get from user_params in .xml
    double sphere_radius = 30.0;
    int nV = 0;

    Cell_Definition* pCD = cell_definitions_by_index[0];
    // Cell_Definition* pCD = cell_definitions_by_name("cancer cell");
    std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
    for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
    {
        std::vector<double> position = {0,0,0};
        // position[0] = Xmin + UniformRandom()*Xrange;
        // position[1] = Ymin + UniformRandom()*Yrange;
        // position[2] = Zmin + UniformRandom()*Zrange;

        // sprinkle some cells near the end of the cylinder/vessel
        double xval = cyl_end_xval + UniformRandom()*sphere_radius;
        double yval = UniformRandom()*sphere_radius;
        double zval = UniformRandom()*sphere_radius;
        position[0] = xval;
        position[1] = yval;
        position[2] = zval;
        std::cout << xval<<", "<< yval<<", "<< zval<<std::endl;

        pC = create_cell( *pCD );
        pC->assign_position( position );

        int m = microenvironment.nearest_voxel_index( position );
        microenvironment(m)[nV] = 32.0;
    }
	std::cout << std::endl;

    //double cyl_x0 = parameters.doubles("cyl_x0"); commented out randy code for axis-aligned vessel
    //double cyl_x1 = parameters.doubles("cyl_x1");

    // added lines to read vessel data from string parameters

		int vessel_number = parameters.strings.parameters.size();
		std::cout <<"vessel_number = "<<vessel_number<<std::endl;

		for( int i = 0 ; i < vessel_number; i++ )
    {
		std::string vessel_n = std::to_string(i+1);
		std::string str_vessel = parameters.strings("vessel_" + vessel_n);
    std::vector<double> vessel;
    std::stringstream s_stream(str_vessel);
    while (s_stream.good()){
        std::string substr;
        std::getline(s_stream,substr,',');
        double value = std::stod(substr);
        vessel.push_back(value);
    }

    double cyl_x0 = vessel[0];
    double cyl_y0 = vessel[1];
    double cyl_z0 = vessel[2];
    double cyl_x1 = vessel[3];
    double cyl_y1 = vessel[4];
    double cyl_z1 = vessel[5];
		double cyl_rad = vessel[6];
    double dist = std::pow((cyl_x1-cyl_x0)*(cyl_x1-cyl_x0) + (cyl_y1-cyl_y0)*(cyl_y1-cyl_y0) + (cyl_z1-cyl_z0)*(cyl_z1-cyl_z0),0.5);
    int nn = dist/10;
    //nn = dist/20;

		// Create new basis vectors given unit vector along cyl
		std::vector<double> n0 = {0,0,0};
		double n0_x = (cyl_x1 - cyl_x0);
		double n0_y = (cyl_y1 - cyl_y0);
		double n0_z = (cyl_z1 - cyl_z0);
		n0[0] = n0_x / dist; // dist is magnitude of n0
		n0[1] = n0_y / dist;
		n0[2] = n0_z / dist;

		// Random point
		std::vector<double> pt_rand_diff = {0,0,0};
		double pt_rand_diff_x = 0.5 - cyl_x0; //srand() - cyl_x0
		double pt_rand_diff_y = 0.6 - cyl_y0; //srand() - cyl_y0
		double pt_rand_diff_z = 0.7 - cyl_z0; //srand() - cyl_z0
		double pt_rand_diff_mag = std::pow(pt_rand_diff_x*pt_rand_diff_x + pt_rand_diff_y*pt_rand_diff_y + pt_rand_diff_z*pt_rand_diff_z,0.5);
		pt_rand_diff[0] = pt_rand_diff_x / pt_rand_diff_mag;
		pt_rand_diff[1] = pt_rand_diff_y / pt_rand_diff_mag;
		pt_rand_diff[2] = pt_rand_diff_z / pt_rand_diff_mag;
		double pt_rand_diff_dot = pt_rand_diff[0] * n0[0] + pt_rand_diff[1] * n0[1] + pt_rand_diff[2] * n0[2];

		std::vector<double> n1 = {0,0,0};
		double n1_x = pt_rand_diff[0] - pt_rand_diff_dot * n0[0];
		double n1_y = pt_rand_diff[1] - pt_rand_diff_dot * n0[1];
		double n1_z = pt_rand_diff[2] - pt_rand_diff_dot * n0[2];
		double n1_mag = std::pow(n1_x*n1_x + n1_y*n1_y + n1_z*n1_z,0.5);
		n1[0] = n1_x / n1_mag;
		n1[1] = n1_y / n1_mag;
		n1[2] = n1_z / n1_mag;

		std::vector<double> n2 = {0,0,0};
		double n2_x = n0_y * n1_z - n0_z * n1_y;
		double n2_y = -(n0_x * n1_z - n0_z * n1_x);
		double n2_z = n0_x * n1_y - n0_y * n1_x;
		double n2_mag = std::pow(n2_x*n2_x + n2_y*n2_y + n2_z*n2_z,0.5);
		n2[0] = n2_x / n2_mag;
		n2[1] = n2_y / n2_mag;
		n2[2] = n2_z / n2_mag;

		double n0n1 = n0[0]*n1[0] + n0[1]*n1[1] + n0[2]*n1[2];
		double n0n2 = n0[0]*n2[0] + n0[1]*n2[1] + n0[2]*n2[2];
		double n1n2 = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2];

		std::cout <<"n0 dot n1 = "<<n0n1<<std::endl;
		std::cout <<"n0 dot n2 = "<<n0n2<<std::endl;
		std::cout <<"n1 dot n2 = "<<n1n2<<std::endl;

		double pi = 3.14159265358979323846;
		double n_cells = 10;
		double del_theta = (2 * pi) / n_cells;

    //cyl_x0 = -200.0;
    pCD = cell_definitions_by_index[1];
    std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
    for( int n = 0 ; n < nn; n++ )
    {
        std::vector<double> center_position = {0,0,0};
        // position[0] = Xmin + UniformRandom()*Xrange;
        // position[1] = Ymin + UniformRandom()*Yrange;
        // position[2] = Zmin + UniformRandom()*Zrange;

        // Place center points along vessel
        double center_x = cyl_x0 + n*((cyl_x1-cyl_x0)/nn);
        double center_y = cyl_y0 + n*((cyl_y1-cyl_y0)/nn);
        double center_z = cyl_z0 + n*((cyl_z1-cyl_z0)/nn);
        center_position[0] = center_x;
        center_position[1] = center_y;
        center_position[2] = center_z;

        std::cout << n << ") " << center_x<<", "<< center_y<<", "<< center_z<<std::endl;

				//Add n cells in a circle around the center point
				for( int i = 0 ; i < n_cells; i++ )
				{
				std::vector<double> position = {0,0,0};
				double xval = center_x + cyl_rad * (cos(i * del_theta) * n1[0] + sin(i * del_theta) * n2[0]);
				double yval = center_y + cyl_rad * (cos(i * del_theta) * n1[1] + sin(i * del_theta) * n2[1]);
				double zval = center_z + cyl_rad * (cos(i * del_theta) * n1[2] + sin(i * del_theta) * n2[2]);
				position[0] = xval;
				position[1] = yval;
				position[2] = zval;

				pC = create_cell( *pCD );
				pC->assign_position( position );
				}
    }
	}


	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();

	return;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; }

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; }
