#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

#define DEBUG_SELDONDATA_INDICES
#define DEBUG_SELDONDATA_IO
#define DEBUG_SELDONDATA_DIMENSION
#define DEBUG_SELDONDATA_MEMORY

#include "SeldonData.hxx"
using namespace SeldonData;

// Permutation of SPACK input.
// 3 parameters are required:
// 'FileSpecies': SPACK input for species names.
// 'FileSpeciesPerm': new SPACK input for species names.
// 'FilePerm': permutation.
int main(int argc, char** argv)
{

  TRY;

  typedef double real;

  int i, j, p, q;
  bool add;
  string temp;

  int Ns;

  ifstream InputSpecies;
  InputSpecies.open(argv[1], ifstream::in);

  getline(InputSpecies, temp);
  getline(InputSpecies, temp);

  InputSpecies >> Ns;

  InputSpecies.close();

  // Output file.
  ofstream Output;
  Output.open(argv[2]);

  string align = "      ";

  // Permutation.

  Grid<int> GridX(Ns);

  Data<int, 1> perm(GridX);

  // Read inputs.

  ifstream Input_perm;
  Input_perm.open(argv[3]);

  for (i = 0; i < Ns; i++)
    Input_perm >> perm(i);

  ifstream Input_Scheme;
  Input_Scheme.open(argv[1]);

  getline(Input_Scheme, temp);
  Output << temp << endl;
  getline(Input_Scheme, temp);
  Output << temp << endl;
  getline(Input_Scheme, temp);
  Output << temp << endl;
  getline(Input_Scheme, temp);
  Output << temp << endl;

  vector<string> Species(Ns);
  for (i = 0; i < Ns; i++)
    getline(Input_Scheme, Species[i]);

  // Fortran file.

  for (i = 0; i < Ns; i++)
    Output << Species[perm(i) - 1] << endl;;

  END;

  return 0;

}
