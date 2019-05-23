%module talos
%{
#include "String.hxx"
#include "Date.hxx"
#include "Files.hxx"
  %}

%include "std_string.i"
%include "std_vector.i"
namespace std
{
  %template(VectorString) vector<string>;
  %template(VectorInt) vector<int>;
  %template(VectorDouble) vector<double>;
}

using namespace std;

%exception
{
  try
    {
      $action
	}
  catch(std::exception& e)
    {
      PyErr_SetString(PyExc_Exception, e.what());
      return NULL;
    }
  catch(std::string& s)
    {
      PyErr_SetString(PyExc_Exception, s.c_str());
      return NULL;
    }
  catch(const char* s)
    {
      PyErr_SetString(PyExc_Exception, s);
      return NULL;
    }
  catch(...)
    {
      PyErr_SetString(PyExc_Exception, "Unknown exception...");
      return NULL;
    }
}

%include "String.hxx"
%include "Date.hxx"
%include "Files.hxx"
