# include "Galaxy.h"

#define N_COMP 6
#define MAXSTARSPERGLOBE 10000
#define MAX(a,b) ((a < b) ?  (b) : (a))
#define MIN(a,b) ((a > b) ?  (b) : (a))


int main (int argc, const char **argv)
{

// Setparameter file and finaloutput filename
	string outfile;
        paramfile params (argv[1],false);
        outfile = params.find<string>("OutFile","demo");

// image related variables

	long numberofparticles=0;
	long maxstarsperglobe;
	float * F_starx;
	float * F_stary;
	float * Red;
	float * Blue;
	float * Green;
	float * III;
	long nx, ny;


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
	printf("======== Processing the Object ========\n");
	printf("=======================================\n");

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

	long nwant,nfinal=0;

	for(int itype=0;itype<N_COMP;itype++)
	  {
	    printf("Generating %s component \n",ComponentsName[itype].c_str());
            long component_type = params.find<long>("Do"+ComponentsName[itype],0);

	    if(component_type == 0)
	      {
		printf("  Component switched off\n");
                nwant = nfinal = npart[itype] = 0;
	      }
	    else
	      {
		nwant = params.find<long>("N"+ComponentsName[itype],0);
		xcomp = new float [nwant];
		ycomp = new float [nwant];
		zcomp = new float [nwant];

		if(component_type == 3 || component_type == 4 || component_type == 5 || component_type == 6 || component_type == 7)
		  {
		    printf("  Reading color & mask images\n");

		    string imagefile_rgb = params.find<string>(ComponentsName[itype]+"FileRGB","NONE");
		    string imagefile_mask = params.find<string>(ComponentsName[itype]+"FileMask","NONE");

		    nx = params.find<long>(ComponentsName[itype]+"xres",1000);
		    ny = params.find<long>(ComponentsName[itype]+"yres",1000);

		    F_starx = new float [nx*ny];
		    F_stary = new float [nx*ny];
		    Red   = new float [nx*ny];
		    Blue  = new float [nx*ny];
		    Green = new float [nx*ny];
		    III   = new float [nx*ny];

		    long infiletype =  params.find<long>("InFileType"+ComponentsName[itype],1);
		    if (infiletype == 0)
		      {
			//numberofstars = ReadBMP(imagefile_rgb, imagefile_maks, 
			//                nx, ny, Rth, Gth, Bth, Red, Green, Blue, III, starx, stary);
		      } 
		    else 
		      {
			numberofparticles = ReadImages(params,imagefile_rgb, imagefile_mask, nx, ny, Red, Green, Blue, III, F_starx, F_stary, nwant);
		      }

		    for(long ii=0;ii<numberofparticles;ii++)
		      {
			xcomp[ii] = F_starx[ii];
			ycomp[ii] = F_stary[ii];
		      }
		  }

		printf("  Generating particle distribution\n");

		switch(component_type)
		  {
		  case 0:
		    printf("    Component switched off\n");
		    break;
		  case 1:
		    printf("    Exponential Spheroid\n");
		    nfinal=GaussRFunc (params, ComponentsName[itype], nwant, xcomp, ycomp, zcomp);
		    break;
		  case 2:
		    printf("    Globular Clusters like\n");
		    maxstarsperglobe = params.find<long>("Nmaxper"+ComponentsName[itype],0);
		    nfinal=GlobularCluster(params, ComponentsName[itype], nwant, maxstarsperglobe, xcomp, ycomp, zcomp);
		    break;
		  case 3:
		    printf("    Generating disk from image\n");
		    nfinal=GaussRDiscFunc (params, ComponentsName[itype], numberofparticles, nwant, xcomp, ycomp, zcomp, nx, ny);
		    break;
		  case 4:
		    printf("    Generating Gas distribution from image\n");
		    nfinal=RDiscFunc (params, ComponentsName[itype], numberofparticles, nwant, xcomp, ycomp, zcomp, III, nx, ny);
		    break;
		  case 5:
		    printf("    Generating spherical symmetric stars distribution from image\n");
		    nfinal=GaussRGlobFunc(params, ComponentsName[itype], numberofparticles, nwant, xcomp, ycomp, zcomp, III, nx, ny);
			//for(int iii=0;iii<nfinal;iii++)cout << xcomp[iii]<< " " << ycomp[iii] << " " << zcomp[iii]  << endl;
		    break;
		  case 6:
		    printf("    Generating distribution from TiRiFiC model\n");
		    nfinal=RDiscFuncTirific (params, ComponentsName[itype], numberofparticles, nwant, xcomp, ycomp, zcomp, III, nx, ny);
		    break;
		  case 7:
		    printf("    Generating distribution from TiRiFiC + DICE models\n");
		    nfinal=RDiscFuncTirificDice (params, ComponentsName[itype], numberofparticles, nwant, xcomp, ycomp, zcomp, III, nx, ny);
                    //for (int ii=1; ii<nwant; ii++)cout << xcomp[ii] << " " << ycomp[ii] << " " << zcomp[ii] << endl;
		    break;
		  case 8:
		    printf("    Generating distribution from DICE Halo model\n");
		    nfinal=RHalo (params, ComponentsName[itype], numberofparticles, nwant, xcomp, ycomp, zcomp, III, nx, ny);
		    break;
		  }

		npart[itype] = nfinal;
		printf("    generated %ld particles instead of %ld\n", nfinal, nwant);

		if(npart[itype] > 0)
		  {
		    vector<float32> color;
		    vector<float32> intensity;
		    vector<float32> hsml;
                    float bnorm = 1.0/float(npart[itype]);

		    cred   = new float [npart[itype]]; 	
		    cgreen = new float [npart[itype]]; 	
		    cblue  = new float [npart[itype]]; 	
		    ciii   = new float [npart[itype]];
                    float brightness_fact = params.find<float>("Brightness"+ComponentsName[itype],1.0);
                    float color_fact      = params.find<float>("ColorScale"+ComponentsName[itype],1.0);
                    float color_red       = params.find<float>("ColorR"+ComponentsName[itype],255.0);
                    float color_green     = params.find<float>("ColorG"+ComponentsName[itype],255.0);
                    float color_blue      = params.find<float>("ColorB"+ComponentsName[itype],255.0);

                    ////float xred = params.find<float>("Red"+ComponentsName[itype],1.0);
                    ////float xgreen = params.find<float>("Green"+ComponentsName[itype],1.0);
                    ////float xblue = params.find<float>("Blue"+ComponentsName[itype],1.0);

		    switch(component_type)
		      {
		      case 0:
			printf("    Nothing to do\n");
			break;
		      case 1:
			printf("    Color read from parameter file (default=white)\n");
/* Claudio's original version
			for (long i=0; i<npart[itype]; i++)
                          {
			    cred[i] = xred;
                            cgreen[i] = xgreen;
                            cblue[i] = xblue;
                            ciii[i] = 1.0;
                          }
*/
			printf("    Assigning brightness (%g)\n", brightness_fact);
			printf("    Assigning RGB colors (%g,%g,%g)\n", color_red, color_green, color_blue);
			for (long i=0; i<npart[itype]; i++)
                          {
			     cred[i] = (color_red + color_fact*bnorm*i);
			     cgreen[i] = (color_green + color_fact*bnorm*i);
			     cblue[i] = (color_blue + color_fact*bnorm*i);
                             //colors vary slightly: factor 0.00001 approx 1/nfinal
                             //most colors should be white = (255,255,255)
                             ciii[i] = brightness_fact;
                          }
			break;
		      case 2:
			printf("    Assigning white color\n");
			for (long i=0; i<npart[itype]; i++)
			  cred[i] = cgreen[i] = cblue[i] = ciii[i] = 1.0;
			break;
		      case 3:
			printf("    Assigning color from image file\n");
			CalculateColours(params, ComponentsName[itype], npart[itype], 
			cred, cgreen, cblue, ciii, Red, Green, Blue, III, xcomp, ycomp, nx, ny);
			break;
		      case 4:
			printf("    Assigning color from image file\n");
			CalculateColours(params, ComponentsName[itype], npart[itype], 
			cred, cgreen, cblue, ciii, Red, Green, Blue, III, xcomp, ycomp, nx, ny);
			//for(int iii=0;iii<100000;iii++)if(ciii[iii] > 0)cout << ciii[iii] << endl;
			break;
		      case 5:
			printf("    Assigning color from image file\n");
			CalculateColours(params, ComponentsName[itype], npart[itype], 
			cred, cgreen, cblue, ciii, Red, Green, Blue, III, xcomp, ycomp, nx, ny);
			break;
		      case 6:
			printf("    Assigning color from image file\n");
			CalculateColours(params, ComponentsName[itype], npart[itype], 
			cred, cgreen, cblue, ciii, Red, Green, Blue, III, xcomp, ycomp, nx, ny);
			break;
                      case 7:
                        printf("    Assigning color from image file\n");
                        CalculateColours(params, ComponentsName[itype], npart[itype],
                        cred, cgreen, cblue, ciii, Red, Green, Blue, III, xcomp, ycomp, nx, ny);
                        break;
		      case 8:
			printf("    Assigning red color\n");
			for (long i=0; i<npart[itype]; i++)
                         {
			  cred[i] = color_red;
                          ciii[i] = brightness_fact;
                          cgreen[i] = cblue[i] = 0.0;
                         }
			break;
		      }
		    //printf("    Rescaling coordinates\n");
                    //float rescale = rescale_coords(params, ComponentsName[itype], xcomp, ycomp, zcomp, npart[itype]);

                    float parsectotirific = params.find<float>("InternToArcsec"+ComponentName,1.0);
		    printf("    Assigning hsml\n");
		    float set_hsml = params.find<float>("hsml"+ComponentsName[itype],0.001); 
                    float hsmlfact = params.find<float>("Adaptivehsml"+ComponentsName[itype],0.0);
                    float hmax = 0.0;
                    float hmin = 1e20;

                    float Imax = 0.0;
                    float Imin = 1e20;
		    for (long i=0; i<npart[itype]; i++)
                      {
                        Imax = MAX(Imax,ciii[i]);
                        Imin = MIN(Imin,ciii[i]);
                      }

                    float I_eff;
		    for (long i=0; i<npart[itype]; i++)
		      {
			xyz.push_back(xcomp[i]);
			xyz.push_back(ycomp[i]);
			xyz.push_back(zcomp[i]);
   
                        if(Imax != Imin)
                        {
                          I_eff = (ciii[i]-Imin)/(Imax-Imin);
                        } else {
                          I_eff = 0.0;
                        }
                        float hsml_eff = parsectotirific*set_hsml/(1.0+hsmlfact*I_eff);
			hsml.push_back(hsml_eff);
                        hmax = MAX(hmax,hsml_eff);
                        hmin = MIN(hmin,hsml_eff);
                        //hmax = MAX(hmax,I_eff);
                        //hmin = MIN(hmin,I_eff);

			intensity.push_back(ciii[i]);

			color.push_back(cred[i]);
			color.push_back(cgreen[i]);
			color.push_back(cblue[i]);
		      }
                    cout << "HMIN, HMAX ------> " << hmin << " " << hmax << endl;

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

		if(component_type == 3 || component_type == 4 || component_type == 6 || component_type == 7)
		  {
		    delete [] F_starx;
		    delete [] F_stary;
		    delete [] Red;
		    delete [] Blue;
		    delete [] Green;
		    delete [] III;
		  }

		delete [] xcomp;
		delete [] ycomp;
		delete [] zcomp;
	      }
	  }


        printf("=====================================\n");
        printf("======= Gas size     : %d\n", npart[0]);
        printf("======= Bulge size   : %d\n", npart[1]);
        printf("======= Halo  size   : %d\n", npart[2]);
        printf("======= Globular size: %d\n", npart[3]);
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
//          arrdata =  H5Dcreate(file_id, field[8].c_str(), H5T_NATIVE_FLOAT, dtotspace, H5P_DEFAULT);
//          H5Dwrite (arrdata, H5T_NATIVE_FLOAT,  H5S_ALL,  H5S_ALL, H5P_DEFAULT, particle_type);
//          H5Dclose(arrdata);
 
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
	  FILE *pFile;
	  pFile = fopen("points.ascii", "w");
          long iaux=0;
	  cout << "WRITING " << xyz.size()/3 << " DATA\n";
          //for(long ii=0; ii<xyz.size()/3; ii=ii+int(xyz.size()/3/100000))
          for(long ii=0; ii<xyz.size()/3; ii=ii+1)
	  {
             iaux=ii*3; 
	     fprintf(pFile, "%f %f %f\n", xyz[iaux],xyz[iaux+1],xyz[iaux+2]);
	  }
	  fclose(pFile);
//#endif

        
}
