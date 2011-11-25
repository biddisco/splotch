# include "Galaxy.h"

#define N_COMP 6
#define MAXSTARSPERGLOBE 10000
#define STAR_FACTOR 10
#define GAS_FACTOR 10

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
	float * F_starx;
	float * F_stary;
	float * zdeep;
        float * particle_type;
	float fibra_range=32.0;
	float * Red;
	float * Blue;
	float * Green;
	float * III;
	long nx, ny;
	float Rth;
	float Gth;
	float Bth;


	printf("=======================================\n");
	printf("===== Start Generating the Galaxy =====\n");
	printf("=======================================\n");

#ifdef HDF5
	cout << "Writing hdf5 data ..." << endl;
#else
	cout << "Writing Gadget format ..." << endl;
	bofstream file(outfile.c_str(),false);
	int32 blocksize=0, blksize=8;
#endif

	printf("=======================================\n");
	printf("======== Processing the Image =========\n");
	printf("=======================================\n");



        nx = params.find<long>("xres",1000);
        ny = params.find<long>("yres",1000);

        Rth = params.find<float>("Rth",0);
        Gth = params.find<float>("Gth",0);
        Bth = params.find<float>("Bth",0);

	float xmax = (float)nx;
	float ymax = (float)nx;
	float zmax = (float)nx;

	F_starx = new float [nx*ny];
	F_stary = new float [nx*ny];
	zdeep = new float [nx*ny];
	Red   = new float [nx*ny];
	Blue  = new float [nx*ny];
	Green = new float [nx*ny];
	III   = new float [nx*ny];

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
	
	   //numberofstars = ReadBMP(params, imagefile, imagefile1, 
           //                nx, ny, Rth, Gth, Bth, Red, Green, Blue, III, starx, stary);
        } else {
           numberofstars = ReadImages(params, nx, ny, Rth, Red, Green, Blue, III, F_starx, F_stary);
        }

// Generate random components

	long totpoints;

	string ComponentsName[N_COMP];
	ComponentsName[0] = "Gas";
        ComponentsName[1] = "Bulge";
        ComponentsName[2] = "Disk";
        ComponentsName[3] = "GCluster";
        ComponentsName[4] = "Stars";
	ComponentsName[5] = "BHs";
	int32 npart[N_COMP];

	vector<float32> xyz;

	float * xcomp;
	float * ycomp;
	float * zcomp;

	float * cred;
	float * cgreen;
	float * cblue;
	float * ciii;

	long nwant,nfinal,counter,ngas,nstars;

	for(int itype=0;itype<N_COMP;itype++)
	  {
	    printf("=====================================\n");
	    switch(itype)
	      {
	      case 0:
		printf("======= Generating Gas disk ==========\n");
		ngas = params.find<long>("DoGas",0);
		if(ngas > 0)
		  nwant = GAS_FACTOR * nx * ny;
		else
		  nwant = 0;
		break;
	      case 1:
		printf("======= Generating the Bulge ========\n");
		nwant = params.find<long>("BulgeSize",0);
		break;
	      case 2:
		printf("======= Generating the Halo =========\n");
		nwant = params.find<long>("DiskSize",0);
		break;
	      case 3:
		printf("======= Generating Globular Clusters \n");
		nwant = params.find<long>("NGlobes",0) * MAXSTARSPERGLOBE;
		break;
	      case 4:
		printf("======= Generating the Stars =========\n");
		nstars=params.find<long>("DoStars",0);
		nwant = nstars * STAR_FACTOR * nx * ny;
		break;
	      case 5:
		printf("======= Nothing to do for BHs ==========\n");
		nwant=0;
		break;
	      }
	    printf("=====================================\n");
	    xcomp = new float [nwant];
	    ycomp = new float [nwant];
	    zcomp = new float [nwant];

	    switch(itype)
	      {
	      case 0:
		if(nwant > 0)
		  {
		    for(long ii=0;ii<numberofstars;ii++)
		      {
			xcomp[ii] = F_starx[ii];
			ycomp[ii] = F_stary[ii];
		      }
		    nfinal=DiscRFunc (params, ComponentsName[0], numberofstars, nwant, xcomp, ycomp, zcomp, xmax, ymax, zmax, zdeep, III, nx, ny);
		  }
		break;
	      case 1:
		nfinal=GaussRFunc (params, ComponentsName[1], nwant, nwant, xcomp, ycomp, zcomp, xmax, ymax, zmax, zdeep, III, nx, ny);
		break;
	      case 2:
		nfinal=GaussRFunc (params, ComponentsName[2], nwant, nwant, xcomp, ycomp, zcomp, xmax, ymax, zmax, zdeep, III, nx, ny);
		break;
	      case 3:
		nfinal=GlobularCluster(params, ComponentsName[3], nwant, MAXSTARSPERGLOBE, xcomp, ycomp, zcomp);
		break;
	      case 4:
		if(nwant > 0)
		  {
		    for(long ii=0;ii<numberofstars;ii++)
		      {
			xcomp[ii] = F_starx[ii];
			ycomp[ii] = F_stary[ii];
		      }
		    nfinal=GaussRFunc (params, ComponentsName[4], numberofstars, nwant, xcomp, ycomp, zcomp, xmax, ymax, zmax, zdeep, III, nx, ny);
		  }
		break;
	      case 5:
		nfinal=0;
		break;
	      }
	    npart[itype] = nfinal;
	    cout << nwant << " -> " << nfinal << endl;
	    if(npart[itype] > 0)
	      {

		vector<float32> color(3*npart[itype]);
		vector<float32> intensity(npart[itype]);
		vector<float32> hsml(npart[itype]);

		cred   = new float [npart[itype]]; 	
		cgreen = new float [npart[itype]]; 	
		cblue  = new float [npart[itype]]; 	
		ciii   = new float [npart[itype]];

		switch(itype)
		  {
		  case 0:
		  case 1:
		  case 2:
		  case 3:
		  case 4:
		  case 5:
		    CalculateColours(npart[itype], cred, cgreen, cblue, ciii, Red, Green, Blue, III, xcomp, ycomp, nx, ny);
		    break;
		  }

		for (long i=counter=0; i<npart[itype]; i++)
		  {
		    xyz.push_back(xcomp[i]);
		    xyz.push_back(ycomp[i]);
		    xyz.push_back(zcomp[i]);

		    hsml[i] = 0.00001;

		    intensity[i] = ciii[i];

		    color[counter++] = cred[i];
		    color[counter++] = cgreen[i];
		    color[counter++] = cblue[i];

		  }

#ifdef HDF5

#else
// write hsml
		string label("HSM"+dataToString(itype));
		file << blksize;
		blocksize = npart[itype]*4 + 8;
		file.put(label.c_str(),4);
		file << blocksize;
		file << blksize;

		file << blocksize-8;
		file.put(&hsml[0],npart[itype]);
		file << blocksize-8;

// write intensity
		label = "INT"+dataToString(itype);
		file << blksize;
		blocksize = npart[itype]*4 + 8;
		file.put(label.c_str(),4);
		file << blocksize;
		file << blksize;

		file << blocksize-8;
		file.put(&intensity[0],npart[itype]);
		file << blocksize-8;

// write color
		label = "COL"+dataToString(itype);
		file << blksize;
		blocksize = 3*npart[itype]*4 + 8;
		file.put(label.c_str(),4);
		file << blocksize;
		file << blksize;

		file << blocksize-8;
		file.put(&color[0],3*npart[itype]);
		file << blocksize-8;
#endif
		delete [] cred;
		delete [] cgreen;
		delete [] cblue;
		delete [] ciii;
	      }

	    delete [] xcomp;
	    delete [] ycomp;
	    delete [] zcomp;

	  }


        printf("=====================================\n");
        printf("======= Gas size     : %d\n", npart[0]);
        printf("======= Bulge size   : %d\n", npart[1]);
        printf("======= Halo  size   : %d\n", npart[2]);
        printf("======= Globiularsize: %d\n", npart[3]);
        printf("======= Stars size   : %d\n", npart[4]);
        printf("======= BHs size     : %d\n", npart[5]);
        printf("=====================================\n");


#ifdef HDF5

#else
// write positions
	string label("POS ");
	file << blksize;
	blocksize = xyz.size()*4 + 8;
	file.put(label.c_str(),4);
	file << blocksize;
	file << blksize;

	file << blocksize-8;
	file.put(&xyz[0],xyz.size());
	file << blocksize-8;
#endif

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

// calculate rho

//        float smooth = params.find<float>("Smooth",0);

     //   CalculateDensity(hsml, rho, xcoord, ycoord, zcoord, nobjects, smooth);


#ifdef HDF5
	cout << "Writing hdf5 data ..." << endl;
// write data in HDF5

	  hid_t file_id = H5Fcreate(outfile.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	  
          hid_t obj_id;
          hid_t dspace;
#define rank 1
          hsize_t dims[rank];
          hsize_t maxdims[rank];
          hid_t dtotspace;
          hid_t arrdata;

	  dims[0] = nobjects;

          //for (long ii=0; ii<nobjects; ii++)floatRed[ii]=(float)ciii[ii];	  

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
          arrdata =  H5Dcreate(file_id, field[5].c_str(), H5T_NATIVE_FLOAT, dtotspace, H5P_DEFAULT);
          H5Dwrite (arrdata, H5T_NATIVE_FLOAT,  H5S_ALL,  H5S_ALL, H5P_DEFAULT, cred);
          H5Dclose(arrdata);
// write green
          arrdata =  H5Dcreate(file_id, field[6].c_str(), H5T_NATIVE_FLOAT, dtotspace, H5P_DEFAULT);
          H5Dwrite (arrdata, H5T_NATIVE_FLOAT,  H5S_ALL,  H5S_ALL, H5P_DEFAULT, cgreen);
          H5Dclose(arrdata);
// write blue
          arrdata =  H5Dcreate(file_id, field[7].c_str(), H5T_NATIVE_FLOAT, dtotspace, H5P_DEFAULT);
          H5Dwrite (arrdata, H5T_NATIVE_FLOAT,  H5S_ALL,  H5S_ALL, H5P_DEFAULT, cblue);
          H5Dclose(arrdata);
// write float particle_type
          arrdata =  H5Dcreate(file_id, field[8].c_str(), H5T_NATIVE_FLOAT, dtotspace, H5P_DEFAULT);
          H5Dwrite (arrdata, H5T_NATIVE_FLOAT,  H5S_ALL,  H5S_ALL, H5P_DEFAULT, particle_type);
          H5Dclose(arrdata);
 
	  H5Fclose(file_id);
#else
// write head
          int32 dummy[64];
	  for(int i=0;i<64;i++)
	    dummy[i]=0;

	  file << blksize;
          blocksize = 256 + 8;
	  file.put("HEAD",4);
	  file << blocksize;
	  file << blksize;

	  file << blocksize-8;
	  file.put(npart,6);
	  file.put(&dummy[0],18);
	  file.put(npart,6);
	  file.put(&dummy[0],64-6-18-6);
	  file << blocksize-8;

	  file.close();
#endif

//#ifdef WRITE_ASCII
	  pFile = fopen("points.ascii", "w");
          long iaux=0;
          for(long ii=0; ii<xyz.size()/3; ii=ii+int(xyz.size()/3/100000))
	  {
             iaux=ii*3; 
	     fprintf(pFile, "%f %f %f\n", xyz[iaux],xyz[iaux+1],xyz[iaux+2]);
	  }
	  fclose(pFile);
//#endif

        
}
