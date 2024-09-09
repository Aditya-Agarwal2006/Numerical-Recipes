#ifndef _NR_TYPES_H_
#define _NR_TYPES_H_

#include <complex>
#include <fstream>
#include "newmat.h"
#include "nrutil.h"
using namespace std;

typedef double DP;

// Vector Types

typedef const NRVec<bool> Vec_I_BOOL;
typedef NRVec<bool> Vec_BOOL, Vec_O_BOOL, Vec_IO_BOOL;

typedef const NRVec<char> Vec_I_CHR;
typedef NRVec<char> Vec_CHR, Vec_O_CHR, Vec_IO_CHR;

typedef const NRVec<unsigned char> Vec_I_UCHR;
typedef NRVec<unsigned char> Vec_UCHR, Vec_O_UCHR, Vec_IO_UCHR;

typedef const NRVec<int> Vec_I_INT;
typedef NRVec<int> Vec_INT, Vec_O_INT, Vec_IO_INT;

typedef const NRVec<unsigned int> Vec_I_UINT;
typedef NRVec<unsigned int> Vec_UINT, Vec_O_UINT, Vec_IO_UINT;

typedef const NRVec<long> Vec_I_LNG;
typedef NRVec<long> Vec_LNG, Vec_O_LNG, Vec_IO_LNG;

typedef const NRVec<unsigned long> Vec_I_ULNG;
typedef NRVec<unsigned long> Vec_ULNG, Vec_O_ULNG, Vec_IO_ULNG;

typedef NRVec<float> Vec_I_SP;
typedef NRVec<float> Vec_SP, Vec_O_SP, Vec_IO_SP;

//typedef const NRVec<DP> Vec_I_DP;
//typedef NRVec<DP> Vec_DP, Vec_O_DP, Vec_IO_DP;

#ifdef use_namespace
typedef const NEWMAT::ColumnVector Vec_I_DP;
typedef NEWMAT::ColumnVector Vec_DP, Vec_O_DP, Vec_IO_DP;
#else
typedef const ColumnVector Vec_I_DP;
typedef ColumnVector Vec_DP, Vec_O_DP, Vec_IO_DP;
#endif

typedef const NRVec<complex<float> > Vec_I_CPLX_SP;
typedef NRVec<complex<float> > Vec_CPLX_SP, Vec_O_CPLX_SP, Vec_IO_CPLX_SP;

typedef const NRVec<complex<DP> > Vec_I_CPLX_DP;
typedef NRVec<complex<DP> > Vec_CPLX_DP, Vec_O_CPLX_DP, Vec_IO_CPLX_DP;

// Matrix Types

typedef const NRMat<bool> Mat_I_BOOL;
typedef NRMat<bool> Mat_BOOL, Mat_O_BOOL, Mat_IO_BOOL;

typedef const NRMat<char> Mat_I_CHR;
typedef NRMat<char> Mat_CHR, Mat_O_CHR, Mat_IO_CHR;

typedef const NRMat<unsigned char> Mat_I_UCHR;
typedef NRMat<unsigned char> Mat_UCHR, Mat_O_UCHR, Mat_IO_UCHR;

typedef const NRMat<int> Mat_I_INT;
typedef NRMat<int> Mat_INT, Mat_O_INT, Mat_IO_INT;

typedef const NRMat<unsigned int> Mat_I_UINT;
typedef NRMat<unsigned int> Mat_UINT, Mat_O_UINT, Mat_IO_UINT;

typedef const NRMat<long> Mat_I_LNG;
typedef NRMat<long> Mat_LNG, Mat_O_LNG, Mat_IO_LNG;

typedef const NRMat<unsigned long> Mat_I_ULNG;
typedef NRMat<unsigned long> Mat_ULNG, Mat_O_ULNG, Mat_IO_ULNG;

typedef const NRMat<float> Mat_I_SP;
typedef NRMat<float> Mat_SP, Mat_O_SP, Mat_IO_SP;

//typedef const NRMat<DP> Mat_I_DP;
//typedef NRMat<DP> Mat_DP, Mat_O_DP, Mat_IO_DP;

#ifdef use_namespace
typedef const NEWMAT::Matrix Mat_I_DP;
typedef NEWMAT::Matrix Mat_DP, Mat_O_DP, Mat_IO_DP;
#else
typedef const Matrix Mat_I_DP;
typedef Matrix Mat_DP, Mat_O_DP, Mat_IO_DP;
#endif

typedef const NRMat<complex<float> > Mat_I_CPLX_SP;
typedef NRMat<complex<float> > Mat_CPLX_SP, Mat_O_CPLX_SP, Mat_IO_CPLX_SP;

typedef const NRMat<complex<DP> > Mat_I_CPLX_DP;
typedef NRMat<complex<DP> > Mat_CPLX_DP, Mat_O_CPLX_DP, Mat_IO_CPLX_DP;

// 3D Matrix Types

typedef const NRMat3d<DP> Mat3D_I_DP;
typedef NRMat3d<DP> Mat3D_DP, Mat3D_O_DP, Mat3D_IO_DP;

// Miscellaneous Types

typedef NRVec<unsigned long *> Vec_ULNG_p;

//typedef NRVec<NRMat<DP> *> Vec_Mat_DP_p;

#ifdef use_namespace
typedef NRVec<NEWMAT::Matrix*> Vec_Mat_DP_p;
#else
typedef NRVec<Matrix*> Vec_Mat_DP_p;
#endif

typedef NRVec<fstream *> Vec_FSTREAM_p;

// additional definition 
#ifdef use_namespace
void print_array(const NEWMAT::ColumnVector &a, const int n_per_row,
   const int col_width);
#else
void print_array(const ColumnVector &a, const int n_per_row,
   const int col_width);
#endif

#endif /* _NR_TYPES_H_ */

// body file: nrtyp_nm.cpp

