#include "reader.h"
#include "ramses_helper_lib.h"
#include <time.h>

//----------------------------------------------------------------------------
// Ramses file reader for particle or amr data (or both)
// Tim Dykes
//
// Use C1, C2, C3, I in parameter file to set the AMR data to visualise. 
// Correct IDs for these can be found in ramses_helper_lib
// 
// Note when using C1+C2+C3 set colour_is_vector0=T, if just one CX variable
// is used then set colour_is_vector0=F.
//
// Particle data drawn using velocity as C1,C2,C3
//
//----------------------------------------------------------------------------


void ramses_reader(paramfile &params, std::vector<particle_sim> &points)
{
	// Check parameter for read mode, points/amr/both
	int mode = params.find<int>("read_mode",0);

	// Get repository, otherwise assume data is in currentdirectory
	std::string repo = params.find<std::string>("infile","");

	// Get all other params
	float active = params.find<float>("active",1);
	float scale3d = params.find<float>("ramses_scale3d",1);

	int ptypes = params.find<int>("ptypes",1);
	std::vector<float> smooth_factor;
	for(int i = 0; i < ptypes; i++)
		smooth_factor.push_back(params.find<float>("smooth_factor"+dataToString(i),1.0));

	std::vector<float> const_intensity;
	for(int i = 0; i < ptypes; i++)
		const_intensity.push_back(params.find<float>("const_intensity"+dataToString(i),-1.0));	

	std::vector<float> intensity;
	for(int i = 0; i < ptypes; i++)
		intensity.push_back(params.find<float>("intensity"+dataToString(i),-1.0));

	std::vector<float> intense_factor;
	for(int i = 0; i < ptypes; i++)
		intense_factor.push_back(params.find<float>("intense_factor"+dataToString(i),1.0));

	std::vector<float> red;
	for(int i = 0; i < ptypes; i++)
		red.push_back(params.find<float>("red"+dataToString(i),-1));

	std::vector<float> green;
	for(int i = 0; i < ptypes; i++)
		green.push_back(params.find<float>("green"+dataToString(i),-1));

	std::vector<float> blue;
	for(int i = 0; i < ptypes; i++)
		blue.push_back(params.find<float>("blue"+dataToString(i),-1));

	srand((unsigned)time(0));

	// Read info file
	info_file info(repo);

	if(mode == 0 || mode == 2)
	{
		// Read amr

		// AMR generated particles are type 0
		int C1 = red[0];
		int C2 = green[0];
		int C3 = blue[0];
		unsigned nlevelmax;

		// Check for constant intensity, if so set i to this value
		bool const_i = false;
		if(const_intensity[0] > -1)
		{
			const_i = true;
			intensity[0] = const_intensity[0];
		}
		// Set default intensity if neither variables are set
		else if(const_intensity[0] == -1 && intensity[0] == -1)
		{
			const_i = true;
			intensity[0] = 1;
		}

		// Check how many colour variables have been set in param file
		// For unset variables, set to density.
		int check = 0;
		if(C1 != -1) ++check;
		else C1 = 0;
		if(C2 != -1) ++check;
		else C2 = 0;
		if(C3 != -1) ++check;
		else C3 = 0;
		if(check == 2) 
		{
			std::cout << "Ramses reader: must set either one or three param file elements red, green, blue, not two. Aborting." << std::endl;
			exit(0);
		}

		if(info.ndim != 3)
		{
			std::cout << "Ramses amr reader: non 3 dimension data detected. Aborting" << std::endl;
			exit(0);
		}


		std::cout << "Reading AMR data... " << std::endl;

		for(unsigned icpu = 0; icpu < info.ncpu; icpu++)
		{
			// Open hydro file and read metadata
			hydro_file hydro(repo, icpu+1);

			// Validate params
			if( (C1 != -1 && (C1 < 0 || C1 >= hydro.meta.nvar)) || (C2 != -1 && (C2 < 0 || C2 >= hydro.meta.nvar)) ||
				(C3 != -1 && (C3 < 0 || C3 >= hydro.meta.nvar)) || (!const_i && (intensity[0] < 0 || intensity[0] >= hydro.meta.nvar)) )
			{
				std::cout << "Trying to read variables that are not in file, there are " << hydro.meta.nvar << " hydro variables in file.";
				std::cout << "Check red,green,blue,intensity parameters. Aborting." << std::endl;
				exit(0);
			}

			// Open amr file + read metadata
			amr_file amr(repo, icpu+1);
			amr.file.SkipRecords(13);
			nlevelmax = amr.meta.nlevelmax; 

			// Read grid numbers
			// Do not account for boundarys. Refer to io_ramses.f90 for boundary support
			F90_Arr2D<unsigned> ngridfile;
			ngridfile.resize(amr.meta.ncpu,amr.meta.nlevelmax);
			amr.file.Read2DArray(ngridfile);

			// Reminder: this line will be different if boundary>0
			amr.file.SkipRecords(3);

			if(info.ordering == "bisection")
				amr.file.SkipRecords(5);
			else 
				amr.file.SkipRecords(4);

			// Storage
			// Not using vectors to avoid performance hit of initialisations on resize
			double* gridvars = 0;
			F90_Arr2D<double> gridcoords;
			F90_Arr3D<double> gridcellhydro;
			F90_Arr2D<unsigned> gridsons;
			int* sonindex = 0;

			// Loop over levels within file
			for(unsigned ilevel = 0; ilevel < amr.meta.nlevelmax; ilevel++)
			{
				// Loop over domains within level (this would be ncpu+nboundary if nboundary>0)
				for(unsigned idomain = 0; idomain < amr.meta.ncpu; idomain++)
				{
					// Check there are grids for this domain
					if(ngridfile(idomain,ilevel) > 0)
					{
						// Allocate arrays
						// Used to tempstore one type of variable for n grids (n grids refers to all grids in this domain and level) (xdp)
						gridvars = new double[ngridfile(idomain,ilevel)];
						// Used to store 3 coords for n grids (xxdp)
						gridcoords.resize(ngridfile(idomain,ilevel), amr.meta.ndim);
						// Used to store set of hydro variables per cell (8) for n grids (vvdp)
						gridcellhydro.resize(ngridfile(idomain,ilevel),amr.meta.twotondim,hydro.meta.nvar);
						// Used to store son indices of n grids (8 indices per grid) (sdp)
						gridsons.resize(ngridfile(idomain,ilevel),amr.meta.twotondim);
						// Tempstore particular cell son index for n grids when reading (idp)
						sonindex = new int[ngridfile(idomain,ilevel)]; 

						// Start reading AMR data
						// Skip grid index, next index and prev index
						amr.file.SkipRecords(3);
						
						// Read centre point of grid for n grids
						for(unsigned idim = 0; idim < amr.meta.ndim; idim++)
						{
							// Read n grid coordinates for 1 dimension and store in gridcoords
							amr.file.Read1DArray(gridvars);
							for(unsigned igrid = 0; igrid < ngridfile(idomain,ilevel); igrid++)
								gridcoords(igrid,idim) = gridvars[igrid];
						}
						// Skip father (1) and neighbour (2ndim) indices
						amr.file.SkipRecords(1+(2*amr.meta.ndim));
						
						// Read son (cell) indices, 8 per grid
						for(unsigned icell = 0; icell < amr.meta.twotondim; icell++)
						{
							amr.file.Read1DArray(sonindex);
							for(unsigned igrid = 0; igrid < ngridfile(idomain,ilevel); igrid++)
								gridsons(igrid,icell) = sonindex[igrid];
						}

						// Skip cpu and refinement maps
						amr.file.SkipRecords(amr.meta.twotondim*2);
					}

					// Start reading hydro data
					hydro.file.SkipRecords(2);

					// Check there are grids for this domain
					if(ngridfile(idomain,ilevel) > 0)
					{
						// Loop over cells (8 per grid)
						for(unsigned icell = 0; icell < amr.meta.twotondim; icell++)
						{
							// Loop over hydro variables
							for(unsigned ivar = 0; ivar < hydro.meta.nvar; ivar++)
							{
								// Read set of hydro variables for a particular cell for n grids
								hydro.file.Read1DArray(gridvars);
								// Loop through n grids for this domain at this level
								for(unsigned igrid = 0; igrid < ngridfile(idomain,ilevel); igrid+=1)
								{
									
									// Store hydro variables in correct location
									gridcellhydro(igrid,icell,ivar) = gridvars[igrid];
									// Conditions: store data only if final loop of hydros (ie, read all hydro data)
									// current domain must match current file..?
									// Current cell must be at highest level of refinement (ie son indexes are 0) or ilevel == lmax
									if((ivar==hydro.meta.nvar-1) && (idomain==icpu) && (gridsons(igrid,icell)==0 || ilevel == amr.meta.nlevelmax))
									{
										double dx = pow(0.5, ilevel);
										int ix, iy, iz;
										iz = icell/4;
										iy = (icell - (4*iz))/2;
										ix = icell - (2*iy) - (4*iz);
										// Calculate absolute coordinates + jitter, and generate particle
										// call ranf(localseed,xx)
                                        // xp(pc,l)=xx*boxlen*dx+xc(l)-boxlen*dx/2
                                        float r = (float)rand()/(float)RAND_MAX;
										particle_sim p;
										p.x = ((((float)rand()/(float)RAND_MAX) * amr.meta.boxlen * dx) +(amr.meta.boxlen * (gridcoords(igrid,0) + (double(ix)-0.5) * dx )) - (amr.meta.boxlen*dx/2)) * scale3d;
										p.y = ((((float)rand()/(float)RAND_MAX) * amr.meta.boxlen * dx) +(amr.meta.boxlen * (gridcoords(igrid,1) + (double(iy)-0.5) * dx )) - (amr.meta.boxlen*dx/2)) * scale3d;
										p.z = ((((float)rand()/(float)RAND_MAX) * amr.meta.boxlen * dx) +(amr.meta.boxlen * (gridcoords(igrid,2) + (double(iz)-0.5) * dx )) - (amr.meta.boxlen*dx/2)) * scale3d;
										// Smoothing length set by box resolution
										// Can be scaled via scale3d if box has been scaled up
										// Can be modified with multiplicative factor smooth_factor
										p.r = (amr.meta.boxlen * scale3d * dx * smooth_factor[0]);
										p.I = (const_i) ? intensity[0] : gridcellhydro(igrid,icell,intensity[0]);
										p.e.r = gridcellhydro(igrid,icell,C1);
										p.e.g = gridcellhydro(igrid,icell,C2);
										p.e.b = gridcellhydro(igrid,icell,C3);
										p.type = 0;
										p.active = active;

										points.push_back(p);
									}
								} // End loop over grids
							} // End loop over hydro vars
						} // End loop over cells
						if(gridvars) delete[] gridvars;
						gridcoords.Delete();
						gridcellhydro.Delete();
						gridsons.Delete();
						if(sonindex) delete[] sonindex;
					} 
				} // End loop over domains
			} // End loop over levels
		} // End loop over files

		std::cout << "Read " << nlevelmax << " levels for " << info.ncpu << " domains." << std::endl;
	}

	// read points
	if(mode == 1 || mode == 2)
	{
		int type = (mode==1) ? 0 : 1;

		int C1 = params.find<int>("red"+dataToString(type),-1);
		int C2 = params.find<int>("green"+dataToString(type),-1);
		int C3 = params.find<int>("blue"+dataToString(type),-1);

		// Check for constant intensity, if so set i to this value
		bool const_i = false;
		if(const_intensity[type] > -1)
		{
			const_i = true;
			intensity[type] = const_intensity[type];
		}
		// Set default intensity if neither variables are set
		else if(const_intensity[type] == -1 && intensity[type] == -1)
		{
			const_i = true;
			intensity[type] = 1;
		}

		// Check how many colour variables have been set in param file
		// For unset variables, set to velocity.
		int check = 0;
		if(C1 != -1) ++check;
		else C1 = 0;
		if(C2 != -1) ++check;
		else C2 = 1;
		if(C3 != -1) ++check;
		else C3 = 2;
		if(check == 2) 
		{
			std::cout << "Ramses particle reader: must set either one or three param file elements red,green,blue, not two. Aborting." << std::endl;
			exit(0);
		}

		// Validate params
		if( (C1 != -1 && (C1 < 0 || C1 > 4)) || (C2 != -1 && (C2 < 0 || C2 > 4)) ||
			(C3 != -1 && (C3 < 0 || C3 > 4)) || (!const_i && (intensity[type] < 0 || intensity[type] > 4)) )
		{
			std::cout << "Trying to read variables that are not readable, check red,green,blue,intensity parameters. Aborting. ";
			exit(0);
		}

		int originalSize = points.size();
		std::cout << "Reading particle data..." << std::endl;
		for(unsigned ifile = 0; ifile < info.ncpu; ifile++)
		{
			// Open file and check header
			part_file part(repo, ifile+1);

			if(part.meta.ndim != 3)
			{
				std::cout << "Ramses particle reader: non 3 dimension data detected. Aborting" << std::endl;
				exit(0);
			}

			// Resize for extra particles
			int previousSize = points.size();
			points.resize(points.size()+part.meta.npart);


			float* fstorage = 0;
			double* dstorage = 0;

			// Check data type, following reads depend on this...
			char dType;
			int dSize;
			part.file.Peek(dSize);
			if(dSize/sizeof(float) == part.meta.npart)
			{
				fstorage = new float[part.meta.npart];
				dType = 'f';
			}
			else if(dSize/sizeof(double) == part.meta.npart)
			{
				dstorage = new double[part.meta.npart];
				dType = 'd';
			}
			else
			{
				std::cout << "Ramses particle reader: Checking data type failed, data size: " << dSize << std::endl;
				std::cout << "npart: " << part.meta.npart << std::endl;
				std::cout << "Aborting." << std::endl;

				exit(0);
			}

			// Read positions
			if(dType=='f')
			{
				part.file.Read1DArray(fstorage);
				for(unsigned i = previousSize; i < points.size(); i++)
					points[i].x = fstorage[i-previousSize] * scale3d;

				part.file.Read1DArray(fstorage);
				for(unsigned i = previousSize; i < points.size(); i++)
					points[i].y = fstorage[i-previousSize] * scale3d;

				part.file.Read1DArray(fstorage);
				for(unsigned i = previousSize; i < points.size(); i++)
					points[i].z = fstorage[i-previousSize] * scale3d;
			}
			else
			{
				part.file.Read1DArray(dstorage);
				for(unsigned i = previousSize; i < points.size(); i++)
					points[i].x = dstorage[i-previousSize] * scale3d;

				part.file.Read1DArray(dstorage);
				for(unsigned i = previousSize; i < points.size(); i++)
					points[i].y = dstorage[i-previousSize] * scale3d;

				part.file.Read1DArray(dstorage);
				for(unsigned i = previousSize; i < points.size(); i++)
					points[i].z = dstorage[i-previousSize] * scale3d;
			}

			// Read appropriate data
			for(int i = 0; i < 5; i++)
			{
				if(C1 == i)
				{
					// If metallicity is requested, skip to metallicity position.
					// Dont request metallicity if it is not in the file...
					if(i==4) part.file.SkipRecords(3);
					if(dType=='f')
					{
						part.file.Read1DArray(fstorage);
						for(unsigned i = previousSize; i < points.size(); i++)
							points[i].e.r = fstorage[i-previousSize];		
					}
					else
					{
						part.file.Read1DArray(dstorage);
						for(unsigned i = previousSize; i < points.size(); i++)
							points[i].e.r = dstorage[i-previousSize];							
					}		
				}
				else if(C2 == i)
				{
					if(i==4) part.file.SkipRecords(3);
					if(dType=='f')
					{
						part.file.Read1DArray(fstorage);
						for(unsigned i = previousSize; i < points.size(); i++)
							points[i].e.g = fstorage[i-previousSize];		
					}
					else
					{
						part.file.Read1DArray(dstorage);
						for(unsigned i = previousSize; i < points.size(); i++)
							points[i].e.g = dstorage[i-previousSize];							
					}	
				}
				else if(C3 == i)
				{
					if(i==4) part.file.SkipRecords(3);
					if(dType=='f')
					{
						part.file.Read1DArray(fstorage);
						for(unsigned i = previousSize; i < points.size(); i++)
							points[i].e.b = fstorage[i-previousSize];		
					}
					else
					{
						part.file.Read1DArray(dstorage);
						for(unsigned i = previousSize; i < points.size(); i++)
							points[i].e.b = dstorage[i-previousSize];							
					}	
				}
				else if(!const_i && intensity[type] == i)
				{
					if(i==4) part.file.SkipRecords(3);
					if(dType=='f')
					{
						part.file.Read1DArray(fstorage);
						for(unsigned i = previousSize; i < points.size(); i++)
							points[i].I = fstorage[i-previousSize] * intense_factor[type];		
					}
					else
					{
						part.file.Read1DArray(dstorage);
						for(unsigned i = previousSize; i < points.size(); i++)
							points[i].I = dstorage[i-previousSize] * intense_factor[type];							
					}	
				}
				// Skip to next record if not final read
				else if(i < 4) part.file.SkipRecord();				
			}


			// Insert param data
			for(unsigned i = previousSize; i < points.size(); i++)
			{
				points[i].r = smooth_factor[type];
				points[i].type = type;
				points[i].active = active;
				if(const_i) points[i].I = intensity[type] * intense_factor[type];
			}

			// Clean up memory
			if(fstorage) delete[] fstorage;
			if(dstorage) delete[] dstorage;
		}

		std::cout << "Read " << points.size() - originalSize << " particles." << std::endl;

	}
	
	if(mode!= 0  && mode != 1 && mode != 2)
	{
		// Explain mode parameter usage and quit
		std::cout << "Reader parameters incorrect. Please set read_mode in parameter file (default 0).\n";
		std::cout << "read_mode=0 for amr data only.\n";
		std::cout << "read_mode=1 for point data only.\n";
		std::cout << "read_mode=2 for both amr + point data.\n";
		std::cout << "Note: You must have amr_ and hydro_ files for opt 0, part_ files for opt1, and all 3 for opt 2." << std::endl;
	}

}