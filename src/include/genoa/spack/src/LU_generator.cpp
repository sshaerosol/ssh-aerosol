#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#define DEBUG_SELDONDATA_INDICES
#define DEBUG_SELDONDATA_IO
#define DEBUG_SELDONDATA_DIMENSION
#define DEBUG_SELDONDATA_MEMORY

#include "SeldonData.hxx"
using namespace SeldonData;

// Returns LU decomposition function and LU solver.
// 6 parameters are required:
// 'FileSpecies': SPACK input for species names.
// 'FileZero': SPACK output (non_zero.dat).
// 'FileLicense': statement for the license.
// 'FileHeaderDec': header of decomposition routine.
// 'FileHeaderSol': header of solving routine.
// 'FileHeaderSolTr': header of solving routine (transposed).
// A 7th parameter can indicate a function suffix to append to the generated
// functions.
int main(int argc, char** argv)
{

  TRY;

  typedef double real;

  // Temporary variables.
  int i, j, p, q;
  bool add;
  string temp, header;

  // Reads the optional function suffix.
  string func_suffix;
  if (argc > 7)
    {
      func_suffix = '_';
      func_suffix += argv[7];
    }

  // Reads the number of species.

  int Ns;

  ifstream InputSpecies;
  InputSpecies.open(argv[1], ifstream::in);

  getline(InputSpecies, temp);
  getline(InputSpecies, temp);

  InputSpecies >> Ns;

  InputSpecies.close();

  // Output files.
  ofstream OutputDec;
  OutputDec.open("LU_decompose.f90", ofstream::out);
  ofstream OutputSol;
  OutputSol.open("LU_solve.f90", ofstream::out);

  // For Fortran...
  string align = "";

  // Matrices M and bool_M.
  RegularGrid<int> GridX(Ns);
  Data<real, 2, int> M(GridX, GridX);
  Data<bool, 2, int> bool_M(GridX, GridX);

  // Reads non zero indices in M.
  // bool_M(i, j) = (M(i, j) !=0).

  for (i = 0; i < Ns; i++)
    for (j = 0; j < Ns; j++)
      if (i != j)
        bool_M(i, j) = false;
      else
        // Diagonal elements cannot be equal to 0.
        bool_M(i, j) = true;

  // Non-zero Extra-diagonal elements.
  ifstream InputNonZero;
  InputNonZero.open(argv[2]);

  // Number of non-zero elements (useless).
  InputNonZero >> i;

  while (InputNonZero.good())
    {
      InputNonZero >> i;
      InputNonZero >> j;
      bool_M(i - 1, j - 1) = true;
    }

  InputNonZero.close();

  // Headers for the license.
  ifstream LicenseHeader;
  LicenseHeader.open(argv[3]);

  string header_license;
  header_license = "";
  while (LicenseHeader.good())
    {
      getline(LicenseHeader, temp);
      header_license += temp + "\n";
    }

  LicenseHeader.close();


  ///////////////////
  // DECOMPOSITION //

  // Headers.
  ifstream DecHeader;
  DecHeader.open(argv[4]);

  header = "";
  while (DecHeader.good())
    {
      getline(DecHeader, temp);
      header += temp + "\n";
    }

  DecHeader.close();

  OutputDec << header_license;

  OutputDec << align << "SUBROUTINE SSH_LU_decompose" << func_suffix << " (ns,M)" \
      << endl;

  OutputDec << endl;
  OutputDec << header;

  OutputDec << align << "IMPLICIT NONE" << endl;

  OutputDec << endl;

  //  OutputDec << "!     -- INCLUDE FILES" << endl;
  //  OutputDec << "!     PARACHEM: parameters for sizes of 'chemical' arrays."
  //   << endl << endl;
  //OutputDec << align << "INCLUDE \'PARACHEM.INC\'" << endl;
  OutputDec << endl;

  OutputDec << align << "INTEGER ns" << endl;
  OutputDec << align << "DOUBLE PRECISION M(ns,ns)" << endl;
  OutputDec << align << "DOUBLE PRECISION temp" << endl;

  OutputDec << endl;
  OutputDec << endl;

  for (i = 0; i < Ns; i++)
    {
      // For L.

      for (p = i; p < Ns; p++)
        {
          add = false;
          for (j = 0; j < i; j++)
            {
              if ((bool_M(j, i))  && (bool_M(p, j)))
                {
                  if (!add)
                    {
                      OutputDec << "!     Lower part." << endl;
                      temp = "temp = M(";
                    }
                  else
                    temp = "temp = temp + M(";

                  // temp += M(p,j) * M(j,i);
                  OutputDec << align + temp
            + to_str(p + 1) + ", " + to_str(j + 1)
            + ") * M("
            + to_str(j + 1) + ", " + to_str(i + 1)
            + ")" << endl;

                  add = true;
                }
            }

          // M(p, i) -= temp;
          if (add)
            {
              OutputDec << align + "M("
        + to_str(p + 1) + ", " + to_str(i + 1)
        + ") = M("
        + to_str(p + 1) + ", " + to_str(i + 1)
        + ") - temp" << endl;
              bool_M(p, i) = true;
            }

        }

      // For U.

      for (q = i + 1; q < Ns; q++)
        {
          add = false;
          for (j = 0; j < i; j++)
            {
              if ((bool_M(i, j)) && (bool_M(j, q)))
                {
                  if (!add)
                    {
                      OutputDec << "!     Upper part." << endl;
                      temp = "temp = M(";
                    }
                  else
                    temp = "temp = temp + M(";

                  // temp += M(i,j) * M(j,q);
                  OutputDec << align + temp
            + to_str(i + 1) + ", " + to_str(j + 1)
            + ") * M("
            + to_str(j + 1) + ", " + to_str(q + 1)
            + ")" << endl;

                  add = true;
                }
            }

          if (add)
            // M(i,q) = ( M(i,q) - temp ) / M(i,i);
            {
              OutputDec << align + "M("
        + to_str(i + 1) + ", " + to_str(q + 1)
        + ") = ( M("
        + to_str(i + 1) + ", " + to_str(q + 1)
        + ") - temp ) / M("
        + to_str(i + 1) + ", " + to_str(i + 1)
        + ")" << endl;

              bool_M(i, q) = true;
            }
          else if (bool_M(i, q))
            // M(i,q) = M(i,q) / M(i,i);
            {
              OutputDec << "!     Upper part." << endl;
              OutputDec << align + "M("
        + to_str(i + 1) + ", " + to_str(q + 1)
        + ") = M("
        + to_str(i + 1) + ", " + to_str(q + 1)
        + ") / M("
        + to_str(i + 1) + ", " + to_str(i + 1)
        + ")" << endl;
            }

        }

      OutputDec << endl;

    }

  OutputDec << endl << align << "END" << endl;

  OutputDec.close();

  // DECOMPOSITION //
  ///////////////////

  ///////////
  // SOLVE //

  ifstream SolHeader;
  SolHeader.open(argv[5]);

  header = "";
  while (SolHeader.good())
    {
      getline(SolHeader, temp);
      header += temp + "\n";
    }

  SolHeader.close();

  OutputSol << header_license;

  OutputSol << align << "SUBROUTINE SSH_LU_solve" << func_suffix << " (ns, M, x)" \
      << endl;

  OutputSol << endl;
  OutputSol << header;

  OutputSol << align << "IMPLICIT NONE" << endl;

  OutputSol << endl;

  // OutputSol << "!     -- INCLUDE FILES" << endl;
  //OutputSol << "!     PARACHEM: parameters for sizes of 'chemical' arrays."
  //<< endl << endl;
  //OutputSol << align << "INCLUDE \'PARACHEM.INC\'" << endl;
  //OutputSol << endl;

  OutputSol << align << "INTEGER ns" << endl;
  OutputSol << align << "DOUBLE PRECISION M(ns,ns)" << endl;
  OutputSol << align << "DOUBLE PRECISION x(ns)" << endl;
  OutputSol << align << "DOUBLE PRECISION temp" << endl;

  OutputSol << endl;
  OutputSol << endl;

  // Forward substitution.

  OutputSol << "!     Forward substitution." << endl;
  OutputSol << endl;

  for (i = 0; i < Ns; i++)
    {
      add = false;
      for (j = 0; j < i; j++)
        if (bool_M(i, j))
          {
            if (!add)
              temp = "temp = M(";
            else
              temp = "temp = temp + M(";

            // temp += M(i, j) * x(j);
            OutputSol << align  + temp
          + to_str(i + 1) + ", " + to_str(j + 1) + ")"
          + " * x(" + to_str(j + 1) + ")"
                      << endl;

            add = true;
          }
      if (add)
        // x(i) = ( x(i) - temp ) / M(i, i);
        OutputSol << align + "x(" + to_str(i + 1)
      + ") = ( x(" + to_str(i + 1) + ") - temp ) / M("
      + to_str(i + 1) + ", " + to_str(i + 1) + ")"
                  << endl;
      else
        // x(i) = x(i) / M(i, i);
        OutputSol << align + "x(" + to_str(i + 1)
      + ") = x(" + to_str(i + 1) + ") / M("
      + to_str(i + 1) + ", " + to_str(i + 1) + ")"
                  << endl;

      OutputSol << endl;
    }

  OutputSol << endl;

  // Backward substitution.

  OutputSol << "!     Backward substitution." << endl;
  OutputSol << endl;

  for (i = Ns - 2; i > -1; i--)
    {
      add = false;
      for (j = i + 1; j < Ns; j++)
        if (bool_M(i, j))
          {
            if (!add)
              temp = "temp = M(";
            else
              temp = "temp = temp + M(";

            // temp += M(i, j) * x(j);
            OutputSol << align  + temp
          + to_str(i + 1) + ", " + to_str(j + 1) + ")"
          + " * x(" + to_str(j + 1) + ")"
                      << endl;

            add = true;
          }
      if (add)
        // x(i) -= temp;
        OutputSol << align + "x(" + to_str(i + 1)
      + ") = x(" + to_str(i + 1) + ") - temp"
                  << endl;

      if (add)
        OutputSol << endl;

    }

  OutputSol << endl << align << "END" << endl;

  OutputSol.close();

  // SOLVE //
  ///////////

  //////////////
  // SOLVE_TR //

    ofstream OutputSolTr;
    OutputSolTr.open("LU_solve_tr.f90", ofstream::out);

    ifstream SolTrHeader;
    SolTrHeader.open(argv[6]);

    header = "";
    while (SolTrHeader.good())
      {
        getline(SolTrHeader, temp);
        header += temp + "\n";
      }

    SolTrHeader.close();

    OutputSolTr << header_license;

    OutputSolTr << align \
        << "SUBROUTINE SSH_LU_solve_tr" << func_suffix << " (ns, M, x)" << endl;

    OutputSolTr << endl;
    OutputSolTr << header;

    OutputSolTr << align << "IMPLICIT NONE" << endl;

    OutputSolTr << endl;

    // OutputSolTr << "!     -- INCLUDE FILES" << endl;
    //OutputSolTr << "!     PARACHEM: parameters for sizes of 'chemical' arrays."
    //<< endl << endl;
    //OutputSolTr << align << "INCLUDE \'PARACHEM.INC\'" << endl;
    //OutputSolTr << endl;

    OutputSolTr << align << "INTEGER ns" << endl;
    OutputSolTr << align << "DOUBLE PRECISION M(ns,ns)" << endl;
    OutputSolTr << align << "DOUBLE PRECISION x(ns)" << endl;
    OutputSolTr << align << "DOUBLE PRECISION temp" << endl;

    OutputSolTr << endl;
    OutputSolTr << endl;

    // Forward substitution.

    OutputSolTr << "!     Forward substitution." << endl;
    OutputSolTr << endl;

    for (i = 1; i < Ns; i++)
      {
        add = false;
        for (j = 0; j < i; j++)
          if (bool_M(j, i))
            {
              if (!add)
                temp = "temp = M(";
              else
                temp = "temp = temp + M(";

              // temp += M(j, i) * x(j);
              OutputSolTr << align  + temp
        + to_str(j + 1) + ", " + to_str(i + 1) + ")"
        + " * x(" + to_str(j + 1) + ")"
                          << endl;

              add = true;
            }
        if (add)
          // x(i) -= temp;
          OutputSolTr << align + "x(" + to_str(i + 1)
        + ") = x(" + to_str(i + 1) + ") - temp"
                      << endl;

        if (add)
          OutputSolTr << endl;
      }

    OutputSolTr << endl;

    // Backward substitution.

    OutputSolTr << "!     Backward substitution." << endl;
    OutputSolTr << endl;

    for (i = Ns - 1; i > -1; i--)
      {
        add = false;
        for (j = i + 1; j < Ns; j++)
          if (bool_M(j, i))
            {
              if (!add)
                temp = "temp = M(";
              else
                temp = "temp = temp + M(";

              // temp += M(j, i) * x(j);
              OutputSolTr << align  + temp
        + to_str(j + 1) + ", " + to_str(i + 1) + ")"
        + " * x(" + to_str(j + 1) + ")"
                          << endl;

              add = true;
            }
        if (add)
          // x(i) = ( x(i) - temp ) / M(i, i);
          OutputSolTr << align + "x(" + to_str(i + 1)
        + ") = ( x(" + to_str(i + 1) + ") - temp ) / M("
        + to_str(i + 1) + ", " + to_str(i + 1) + ")"
                      << endl;
        else
          // x(i) = x(i) / M(i, i);
          OutputSolTr << align + "x(" + to_str(i + 1)
        + ") = x(" + to_str(i + 1) + ") / M("
        + to_str(i + 1) + ", " + to_str(i + 1) + ")"
                      << endl;

        OutputSolTr << endl;

      }

    OutputSolTr << endl << align << "END" << endl;;

    OutputSolTr.close();

  // SOLVE_TR //
  //////////////

  END;

  return 0;

}
