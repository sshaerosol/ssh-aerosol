#include <iostream>
#include <fstream>

using namespace std;

#define DEBUG_SELDONDATA_INDICES
#define DEBUG_SELDONDATA_IO
#define DEBUG_SELDONDATA_DIMENSION
#define DEBUG_SELDONDATA_MEMORY

#include "SeldonData.hxx"
using namespace SeldonData;

// Returns the permutation in file 'FilePerm'.
// 3 parameters are required:
// 'FileSpecies': SPACK input for species names.
// 'FileZero': SPACK output (non_zero.dat).
// 'FilePerm': where to store the permutation.
int main(int argc, char** argv)
{

  TRY;

  typedef double real;

  int i, j, k, p, q;
  bool add;
  string temp;

  // 'main' parameters.
  string species = argv[1];
  string non_zero = argv[2];

  // Reads the number of species.

  int Ns;

  ifstream InputSpecies;
  InputSpecies.open(species.c_str(), ifstream::in);

  getline(InputSpecies, temp);
  getline(InputSpecies, temp);

  InputSpecies >> Ns;

  InputSpecies.close();

  // Reads non zero indices in M.
  // bool_M(i, j) = (M(i, j)!=0).
  Grid<int> GridX(Ns);
  Data<bool, 1, int> C(GridX), R(GridX);
  Data<bool, 2, int> D(GridX, GridX);

  Data<int, 1> perm(GridX);

  D.Fill(false);
  // Diagonal elements cannot be equal to 0.
  for (i = 0; i < Ns; i++)
    D(i, i) = true;

  // Non-zero Extra-diagonal elements.
  ifstream InputNonZero;
  InputNonZero.open(non_zero.c_str(), ifstream::in);

  // Number of non-zero elements (useless).
  InputNonZero >> i;

  while (InputNonZero.good())
    {
      InputNonZero >> i;
      InputNonZero >> j;
      D(i - 1, j - 1) = true;
    }

  InputNonZero.close();

  // Finds the permutation.

  for (i = 0; i < Ns; i++)
    perm(i) = i + 1;

  for (k = 0; k < Ns; k++)
    {
      Grid<int> SubGrid(Ns - k);
      Data<bool, 2, int> P(SubGrid, SubGrid);

      P.SubData(D, Range(k, Ns - 1), Range(k, Ns - 1));

      int sum_min = Ns * Ns;
      int i_min = 0;
      for (i = 0; i < Ns - k; i++)
        {
          int sum = -1;
          for (j = 0; j < Ns - k; j++)
            if (P(i, j)) sum++;
          if (sum_min > sum * sum)
            {
              sum_min = sum * sum;
              i_min = i;
            }
        }

      i = k + i_min;

      // Permutation i <--> k.

      C.SubData(D, Range::all(), i);
      for (j = 0; j < Ns; j++)
        {
          D(j, i) = D(j, k);
          D(j, k) = C(j);
        }

      R.SubData(D, i, Range::all());
      for (j = 0; j < Ns; j++)
        {
          D(i, j) = D(k, j);
          D(k, j) = R(j);
        }

      j = perm(i);
      perm(i) = perm(k);
      perm(k) = j;

      i = k;

      // Simulating LU decomposition.

      for (p = i; p < Ns; p++)
        {
          add = false;
          for (j = 0; j < i; j++)
            add = add || (D(p, j) && D(j, i));
          if (add)
            D(p, i) = true;
        }

      for (q = i + 1; q < Ns; q++)
        {
          add = false;
          for (j = 0; j < i; j++)
            add = add || (D(i, j) && D(j, q));
          if (add)
            D(i, q) = true;
        }

    }

  // Stores the permutation.

  ofstream Output;
  Output.open(argv[3]);

  for (i = 0; i < Ns; i++)
    Output << perm(i) << endl;

  Output.close();

  END;

  return 0;

}
