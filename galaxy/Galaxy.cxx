# include "Galaxy.h"

int main (int argc, const char **argv)
{

        FILE * pFile;

// Setparameter file and finaloutput filename
	string outfile;
        paramfile params (argv[1],false);
        outfile = params.find<string>("OutFile","demo");

// image related variables

	string imagefile1;
	string imagefile;
	long numberofstars;
	float * starx;
	float * stary;
	float * starz;
	float * F_starx;
	float * F_stary;
	float * zdeep;
        float * particle_type;
	float fibra_range=32.0;
	unsigned int * Red;
	unsigned int * Blue;
	unsigned int * Green;
	unsigned int * III;
	unsigned int * F_Red;
	unsigned int * F_Blue;
	unsigned int * F_Green;
	long nx, ny;
	unsigned int Rth;
	unsigned int Gth;
	unsigned int Bth;
	long star_factor = 10;


	printf("=======================================\n");
	printf("===== Start Generating the Galaxy =====\n");
	printf("=======================================\n");
	printf("=======================================\n");
	printf("======== Processing the Image =========\n");
	printf("=======================================\n");



        nx = params.find<long>("xres",1000);
        ny = params.find<long>("yres",1000);

        Rth = params.find<uint>("Rth",0);
        Gth = params.find<uint>("Gth",0);
        Bth = params.find<uint>("Bth",0);

	float xmax = (float)nx;
	float ymax = (float)nx;
	float zmax = (float)nx;

	starx = new float [star_factor*nx*ny];
	stary = new float [star_factor*nx*ny];
	starz = new float [star_factor*nx*ny];
	F_starx = new float [nx*ny];
	F_stary = new float [nx*ny];
	zdeep = new float [nx*ny];
	Red   = new unsigned int [nx*ny];
	Blue  = new unsigned int [nx*ny];
	Green = new unsigned int [nx*ny];
	III   = new unsigned int [nx*ny];
	F_Red   = new unsigned int [nx*ny];
	F_Blue  = new unsigned int [nx*ny];
	F_Green = new unsigned int [nx*ny];

        int infiletype = params.find<int>("InFileType",0);

        if (infiletype == 0)
        {
// BE CAREFUL: I choose nx as basic normalization factor for distances
// FURTHERMORE I assume that galaxy is ALWAYS centered on (0,0,0)

           string galaxyR = params.find<string>("GalaxyFileR","NONE");
           string galaxyG = params.find<string>("GalaxyFileG","NONE");
           string galaxyB = params.find<string>("GalaxyFiler","NONE");
           string galaxyI = params.find<string>("GalaxyFileI","NONE");
           imagefile = galaxyR;
           imagefile1 = galaxyI;
	
	   numberofstars = ReadBMP(params, imagefile, imagefile1, 
                           nx, ny, Rth, Gth, Bth, Red, Green, Blue, III, starx, stary);
        } else {
           numberofstars = ReadImages(params, nx, ny, Rth, Red, Green, Blue, starx, stary);
        }

// Generate random components

	long totpoints;
	const int N_COMP = 4;
	string ComponentsName[N_COMP];
        ComponentsName[0] = "Bulge";
        ComponentsName[1] = "Disk";
        ComponentsName[2] = "Arms";
        ComponentsName[3] = "GCluster";

	printf("=====================================\n");
	printf("======= Generating the Bulge ========\n");
	printf("=====================================\n");

	long bulgesize;
	float * bulgex;
	float * bulgey;
	float * bulgez;
	long nbulge;

        bulgesize = params.find<long>("BulgeSize",0);

	bulgex = new float [bulgesize];
	bulgey = new float [bulgesize];
	bulgez = new float [bulgesize];
	totpoints = bulgesize;

        nbulge=GaussRFunc (params, ComponentsName[0], bulgesize, totpoints, bulgex, bulgey, bulgez, xmax, ymax, zmax, zdeep, III, nx, ny);

        printf("=====================================\n");
        printf("======= Generating the Halo =========\n");
        printf("=====================================\n");

        long halosize;
        float * halox;
        float * haloy;
        float * haloz;
        long nhalo;

        halosize = params.find<long>("DiskSize",0);

        halox = new float [halosize];
        haloy = new float [halosize];
        haloz = new float [halosize];
        totpoints = halosize;

        nhalo=GaussRFunc (params, ComponentsName[1], halosize, totpoints, halox, haloy, haloz, xmax, ymax, zmax, zdeep, III, nx, ny);

        printf("=====================================\n");
        printf("======= Generating Globular Clusters \n");
        printf("=====================================\n");

        float * globex;
        float * globey;
        float * globez;
        long nglobes;
	long maxstarsperglobe=10000;
	long globesize;
	long globenumber;

        nglobes = params.find<long>("NGlobes",0);
	globesize = nglobes*maxstarsperglobe;

        globex = new float [globesize];
        globey = new float [globesize];
        globez = new float [globesize];
        totpoints = globesize;

        globenumber=GlobularCluster(params, ComponentsName[3], nglobes, maxstarsperglobe, globex, globey, globez);

	float * gx;
	float * gy;
	float * gz;
	gx = new float [globenumber];
	gy = new float [globenumber];
	gz = new float [globenumber];
	for (long kk=0; kk<globenumber; kk++)
	{
	   gx[kk] = globex[kk];
	   gy[kk] = globey[kk];
	   gz[kk] = globez[kk];
	}

	delete [] globex;
	delete [] globey;
	delete [] globez;


        printf("=====================================\n");
        printf("======= Generating the Arms =========\n");
        printf("=====================================\n");

	long nstars;
        totpoints = star_factor*nx*ny;

        nstars=GaussRFunc (params, ComponentsName[2], numberofstars, totpoints, starx, stary, starz, xmax, ymax, zmax, zdeep, III, nx, ny);

// BUILD UP THE COMPLETE DATASET

	long nobjects=nhalo+nbulge+nstars+globenumber;

        printf("=====================================\n");
        printf("======= Bulge size   : %d\n", nbulge);
        printf("======= Halo  size   : %d\n", nhalo );
        printf("======= Stars size   : %d\n", nstars);
        printf("======= Globular size: %d\n", globenumber);
        printf("=====================================\n");

	string field[NUM_OF_FIELDS];
        field[0] = "Xpos";
        field[1] = "Ypos";
        field[2] = "Zpos";
        field[3] = "Rho";
        field[4] = "HSML";
        field[5] = "Red";
        field[6] = "Green";
        field[7] = "Blue";
        field[8] = "Type";
        field[9] = "floatGreen";
        field[10] = "floatBlue";

        float * rho;
        float * hsml;
	float * xcoord;
	float * ycoord;
	float * zcoord;
	float * floatRed;
        unsigned int * cred;
        unsigned int * cgreen;
        unsigned int * cblue;
        unsigned int * ciii;

	xcoord = new float [nobjects]; 	
	ycoord = new float [nobjects]; 	
	zcoord = new float [nobjects]; 	
	rho    = new float [nobjects]; 	
	hsml   = new float [nobjects]; 	
	floatRed   = new float [nobjects]; 	
	particle_type  = new float [nobjects]; 	

	long totcounter = 0;

	for (long i=0; i<nstars; i++)
	{
	   xcoord[totcounter] = starx[i];
	   ycoord[totcounter] = stary[i];
	   zcoord[totcounter] = starz[i];
	   rho[totcounter]    = 0.0;
	   hsml[totcounter]   = 0.0;
           particle_type[totcounter] = 0.0;
	   totcounter++;
	}
	delete [] starx;
	delete [] stary;
	delete [] starz;

        for (long i=0; i<nbulge; i++)
        {
           xcoord[totcounter] = bulgex[i];
           ycoord[totcounter] = bulgey[i];
           zcoord[totcounter] = bulgez[i];
           rho[totcounter]    = 0.0;
           hsml[totcounter]   = 0.0;
           particle_type[totcounter] = 1.0;
	   totcounter++;
        }
        delete [] bulgex;
        delete [] bulgey;
        delete [] bulgez;

        for (long i=0; i<nhalo; i++)
        {
           xcoord[totcounter] = halox[i];
           ycoord[totcounter] = haloy[i];
           zcoord[totcounter] = haloz[i];
           rho[totcounter]    = 0.0;
           hsml[totcounter]   = 0.0;
           particle_type[totcounter] = 2.0;
           totcounter++;
        }
        delete [] halox;
        delete [] haloy;
        delete [] haloz;

        for (long i=0; i<globenumber; i++)
        {
           xcoord[totcounter] = gx[i];
           ycoord[totcounter] = gy[i];
           zcoord[totcounter] = gz[i];
           rho[totcounter]    = 0.0;
           hsml[totcounter]   = 0.0;
           particle_type[totcounter] = 3.0;
           totcounter++;
        }
        delete [] gx;
        delete [] gy;
        delete [] gz;


	printf("=====================================\n");
	printf("======== Calculating Density ========\n");
	printf("=====================================\n");


// calculate rho

        float smooth = params.find<float>("Smooth",0);

     //   CalculateDensity(hsml, rho, xcoord, ycoord, zcoord, nobjects, smooth);

	printf("=====================================\n");
	printf("======== Calculating Colours ========\n");
	printf("=====================================\n");

	cred   = new unsigned int [nobjects]; 	
	cgreen = new unsigned int [nobjects]; 	
	cblue  = new unsigned int [nobjects]; 	
	ciii   = new unsigned int [nobjects]; 	

	CalculateColours(nobjects, cred, cgreen, cblue, ciii, Red, Green, Blue, III, xcoord, ycoord, nx, ny);

// write data in HDF5

	  hid_t file_id = H5Fcreate(outfile.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	  
          hid_t obj_id;
          hid_t dspace;
          int rank = 1;
          hsize_t dims[rank];
          hsize_t maxdims[rank];
          hid_t dtotspace;
          hid_t arrdata;

	  dims[0] = nobjects;

          for (long ii=0; ii<nobjects; ii++)floatRed[ii]=(float)ciii[ii];	  

// create geometry
          dtotspace = H5Screate_simple (rank, dims, NULL);
// write x coords
          arrdata =  H5Dcreate(file_id, field[0].c_str(), H5T_NATIVE_FLOAT, dtotspace, H5P_DEFAULT);
          H5Dwrite (arrdata, H5T_NATIVE_FLOAT,  H5S_ALL,  H5S_ALL, H5P_DEFAULT, xcoord);
          H5Dclose(arrdata);
// write y coords
          arrdata =  H5Dcreate(file_id, field[1].c_str(), H5T_NATIVE_FLOAT, dtotspace, H5P_DEFAULT);
          H5Dwrite (arrdata, H5T_NATIVE_FLOAT,  H5S_ALL,  H5S_ALL, H5P_DEFAULT, ycoord);
          H5Dclose(arrdata);
// write z coords
          arrdata =  H5Dcreate(file_id, field[2].c_str(), H5T_NATIVE_FLOAT, dtotspace, H5P_DEFAULT);
          H5Dwrite (arrdata, H5T_NATIVE_FLOAT,  H5S_ALL,  H5S_ALL, H5P_DEFAULT, zcoord);
          H5Dclose(arrdata);
/*
// write rho
          arrdata =  H5Dcreate(file_id, field[3].c_str(), H5T_NATIVE_FLOAT, dtotspace, H5P_DEFAULT);
          H5Dwrite (arrdata, H5T_NATIVE_FLOAT,  H5S_ALL,  H5S_ALL, H5P_DEFAULT, rho);
          H5Dclose(arrdata);
// write hsml
          arrdata =  H5Dcreate(file_id, field[4].c_str(), H5T_NATIVE_FLOAT, dtotspace, H5P_DEFAULT);
          H5Dwrite (arrdata, H5T_NATIVE_FLOAT,  H5S_ALL,  H5S_ALL, H5P_DEFAULT, hsml);
          H5Dclose(arrdata);
*/
// write red
          arrdata =  H5Dcreate(file_id, field[5].c_str(), H5T_NATIVE_UINT, dtotspace, H5P_DEFAULT);
          H5Dwrite (arrdata, H5T_NATIVE_UINT,  H5S_ALL,  H5S_ALL, H5P_DEFAULT, cred);
          H5Dclose(arrdata);
// write green
          arrdata =  H5Dcreate(file_id, field[6].c_str(), H5T_NATIVE_UINT, dtotspace, H5P_DEFAULT);
          H5Dwrite (arrdata, H5T_NATIVE_UINT,  H5S_ALL,  H5S_ALL, H5P_DEFAULT, cgreen);
          H5Dclose(arrdata);
// write blue
          arrdata =  H5Dcreate(file_id, field[7].c_str(), H5T_NATIVE_UINT, dtotspace, H5P_DEFAULT);
          H5Dwrite (arrdata, H5T_NATIVE_UINT,  H5S_ALL,  H5S_ALL, H5P_DEFAULT, cblue);
          H5Dclose(arrdata);
// write float particle_type
          arrdata =  H5Dcreate(file_id, field[8].c_str(), H5T_NATIVE_FLOAT, dtotspace, H5P_DEFAULT);
          H5Dwrite (arrdata, H5T_NATIVE_FLOAT,  H5S_ALL,  H5S_ALL, H5P_DEFAULT, particle_type);
          H5Dclose(arrdata);
 
	  H5Fclose(file_id);


	  pFile = fopen("points.ascii", "w");
          for(long ii=0; ii<nobjects; ii=ii+int(nobjects/10000))
	  {
	     fprintf(pFile, "%f %f %f %f %f\n", xcoord[ii],ycoord[ii],zcoord[ii],floatRed[ii],particle_type[ii]);
	  }
	  fclose(pFile);

	  delete [] xcoord;
	  delete [] ycoord;
	  delete [] zcoord;

        
}
