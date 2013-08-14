#ifndef F90_UNFORMATTEDIO_H
#define F90_UNFORMATTEDIO_H

//----------------------------------------------------------------------------
// Helper library for I/O with Fortran unformatted write files
// Tim Dykes
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Simple templates for 2/3D arrays used for reading file
//----------------------------------------------------------------------------
template<typename T>
class F90_Arr2D{
public: 
	F90_Arr2D() { arr = NULL;}
	void resize(unsigned x, unsigned y)
	{
		xdim = x;
		ydim = y;
		arr = new T[x*y];
	}

	void Delete(){ if(arr != NULL) delete[] arr;}

	T& operator ()(unsigned x, unsigned y) {return arr[x*ydim+y]; }

	unsigned xdim;
	unsigned ydim;

private: 
	T* arr;
};

template<typename T>
class F90_Arr3D{
public: 
	F90_Arr3D() { arr = NULL;}
	void resize(unsigned x, unsigned y, unsigned z)
	{
		xdim = x;
		ydim = y;
		zdim = z;
		arr = new T[x*y*z];
	}

	void Delete(){ if(arr != NULL) delete[] arr;}

	T& operator ()(unsigned x, unsigned y, unsigned z) {return arr[(x*ydim+y)*zdim+z]; }

	unsigned xdim;
	unsigned ydim;
	unsigned zdim;

private: 
	T* arr;
};


//----------------------------------------------------------------------------
// Unformatted file class, allows to read scalar and unknown size 1D arrays
// 2D and 3D arrays must be of a known size
// Allows skipping n records of any size (equivalent to READ(unit) ) 
//----------------------------------------------------------------------------

// Size of prefix and suffix record delimiters in bytes
// These delimiters indicate bytelength of record
#define PREPOST_DATA 4


class F90_UnformattedIO {
public:
	F90_UnformattedIO() {}

	~F90_UnformattedIO() {	file.close(); }

	// Open fortran unformatted file
	void Open(std::string filename) 
	{
		file.open(filename.c_str(), std::ios::in);
		if(!file)
			std::cout << "Failed to open fortran file: " << filename << std::endl;
	}

	template <typename T> 
	void ReadScalar(T& scalar)
	{
		// Read prefix, data and suffix. Check prefix and suffix = sizeof data expected
		unsigned pre, post;

		file.read((char*)&pre, PREPOST_DATA);
		file.read((char*)&scalar,sizeof(T));
		file.read((char*)&post,PREPOST_DATA);

		if((pre != sizeof(T)) || (pre != post)) {
			std::cout << "Failed scaler read from fortran file..." << std::endl;
		}
	}

	// Read 1d array of unknown size
	template <typename T>
	void Read1DArray(T* arr)
	{
		unsigned pre, post;
		file.read((char*)&pre, PREPOST_DATA);

		for(unsigned n = 0; n < (pre/sizeof(T)); n++)
			file.read((char*)&arr[n], sizeof(T));

		file.read((char*)&post, PREPOST_DATA);
		if(pre!=post)
			std::cout << "Failed read fortran 1d array."<< std::endl;
	}

	// Read 2d array of predefined size arr[xdim][ydim]
	// Note: F90_Arr2D/3D must be resized to the expected read size before reading.
	// This is due to no indication within record as to multidimensional array length
	template <typename T>
	void Read2DArray(F90_Arr2D<T> arr)
	{
		unsigned pre, post;
		file.read((char*)&pre, PREPOST_DATA);
		
		for(unsigned n = 0; n < arr.ydim; n++)
			for(unsigned m = 0; m < arr.xdim; m++) 
				file.read((char*)&arr(m,n), sizeof(T));

		file.read((char*)&post, PREPOST_DATA);
		if(pre!=post)
			std::cout << "Failed read fortran 2d array..." << std::endl;
	}

	// Read 3d array of predefined size arr[xdim][ydim][zdim]
	template <typename T>
	void Read3DArray(F90_Arr3D<T> arr)
	{
		unsigned pre, post;
		file.read((char*)&pre, PREPOST_DATA);

		for(unsigned l = 0; l < arr.zdim; l++)
			for(unsigned n = 0; n < arr.ydim; n++)
				for(unsigned m = 0; m < arr.xdim; m++) {
					file.read((char*)&arr(m,n,l), sizeof(T));
				}

		file.read((char*)&post, PREPOST_DATA);

		if(pre!=post)
			std::cout << "Failed read fortran 3d array..." << std::endl;
	}

	// Skip single record (scalar or array)
	void SkipRecord()
	{
		unsigned pre, post;
		try{
			file.read((char*)&pre, PREPOST_DATA);
			file.seekg(pre,std::ios::cur);
			file.read((char*)&post,PREPOST_DATA);
		}catch(...){
			throw std::runtime_error("Error seeking end of record");
		}

		if(pre != post)
			std::cout << "Failed record skip" << std::endl;
	}

	// Skip N records (scalar or array)
	void SkipRecords(unsigned n)
	{
		for(unsigned i = 0; i < n; i++)
		{
			unsigned pre, post;
			try{
				file.read((char*)&pre, PREPOST_DATA);
				file.seekg(pre,std::ios::cur);
				file.read((char*)&post,PREPOST_DATA);
			}catch(...){
				throw std::runtime_error("Error seeking end of record");
			}

			if(pre != post)
				std::cout << "Failed record skip" << std::endl;
		}
	}

	// Peek at next data item without moving on
	// Can be used to check delimiter contents
	template <typename T>
	void Peek(T& data)
	{
		file.read((char*)&data, sizeof(T));
		file.seekg(-sizeof(T),std::ios::cur);
	}

private:
	std::ifstream file;

};


#endif