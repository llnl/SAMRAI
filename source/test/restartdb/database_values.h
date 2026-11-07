/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2026 Lawrence Livermore National Security, LLC
 * Description:   Global values for the restart tests
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/DatabaseBox.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"

#include <string>

using namespace SAMRAI;

// Number of (non-abortive) failures.
int number_of_failures = 0;

float db_float_val = (float)3.14159;
int db_int_val = 4;

dcomplex arraydb_dcomplexArray0 = dcomplex(1, 2);
dcomplex arraydb_dcomplexArray1 = dcomplex(2, 3);
dcomplex arraydb_dcomplexArray2 = dcomplex(3, 4);
bool arraydb_boolArray0 = true;
bool arraydb_boolArray1 = false;
bool arraydb_boolArray2 = false;
int arraydb_intArray0 = 0;
int arraydb_intArray1 = 1;
int arraydb_intArray2 = 2;
int arraydb_intArray3 = 3;
int arraydb_intArray4 = 4;
std::string arraydb_stringArray0 = "This is 1 test.";
std::string arraydb_stringArray1 = "This is 2nd test.";
std::string arraydb_stringArray2 = "This is a long 3rd test.";
float arraydb_floatArray0 = 0 * 1.2;
float arraydb_floatArray1 = (float)(1 * 1.2);
float arraydb_floatArray2 = (float)(2 * 1.2);
float arraydb_floatArray3 = (float)(3 * 1.2);
float arraydb_floatArray4 = (float)(4 * 1.2);
double arraydb_doubleArray0 = 0 * 1.111111;
double arraydb_doubleArray1 = 1 * 1.111111;
double arraydb_doubleArray2 = 2 * 1.111111;
double arraydb_doubleArray3 = 3 * 1.111111;
double arraydb_doubleArray4 = 4 * 1.111111;
double arraydb_doubleArray5 = 5 * 1.111111;
char arraydb_charArray0 = 'a';
char arraydb_charArray1 = 'b';
tbox::DatabaseBox arraydb_boxArray0;
tbox::DatabaseBox arraydb_boxArray1;
tbox::DatabaseBox arraydb_boxArray2;

float scalardb_float1 = (float)1.111;
float scalardb_float2 = (float)2.222;
float scalardb_float3 = (float)3.333;
double scalardb_full_thisDouble = 123.456;
dcomplex scalardb_full_thisComplex = dcomplex(2.3, 4.5);
int scalardb_full_thisInt = 89;
float scalardb_full_thisFloat = (float)9.9;
bool scalardb_full_thisBool = true;
std::string scalardb_full_thisString = "This is a test.";
char scalardb_full_thisChar = 'q';
int ilo[3] = { 0, 0, 0 };
int ihi[3] = { 1, 1, 1 };
tbox::DatabaseBox scalardb_full_thisBox(tbox::Dimension(3), ilo, ihi);

hier::IntVector intVector0(tbox::Dimension(2), 0);
hier::IntVector intVector1(tbox::Dimension(2), 1);
hier::IntVector intVector2(tbox::Dimension(2), 1);
