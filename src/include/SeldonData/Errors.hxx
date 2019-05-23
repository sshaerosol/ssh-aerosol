// Copyright (C) 2003-2007, INRIA
// Author(s): Vivien Mallet
//
// This file is part of SeldonData library, used for data processing.
//
// SeldonData is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// SeldonData is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// For more information, visit the SeldonData home page:
//      http://vivienmallet.net/lib/seldondata/

#ifndef FILE_SELDONDATA_ERRORS_HXX

#include <iostream>
#include <sstream>
#include <string>

namespace SeldonData
{

  ///////////
  // ERROR //
  ///////////

  //! Base class.
  class Error : public std::exception
  {

  protected:
    //! Formatted message for the end user.
    string message_;

  public:
    //! Constructor.
    Error (const string& function = "", const string& comment = "",
           const string& description = "An unknown error occurred")
    {
      cerr << "ERROR!" << endl;

      stringstream s;
      s << description;
      if (!function.empty ())
        s << " in '" << function << "'";
      s << "." << endl;
      if (!comment.empty ())
        s << "   " << comment << endl;
      s << endl;

      message_ = s.str ();

      cerr << message_ << endl;
      cerr.flush();
    }

    //! Destructor.
    virtual
    ~Error () throw ()
    {
    }

    //! Displays error description.
    void
    What ()
    {
      cerr << message_;
    }

    //! Returns an error description, overridden from 'std::exception'.
    virtual const char*
    what ()
    {
      return message_.c_str ();
    }

  };

  //////////////
  // NOMEMORY //
  //////////////

  //! No memory available.
  class NoMemory : public Error
  {

  public:
    //! Constructor.
    NoMemory (const string& function = "", const string& comment = "") :
        Error (function, comment, "No more memory is available")
    {
      // Nothing
    }

  };

  //////////////
  // WRONGDIM //
  //////////////

  //! Wrong dimension.
  /*!
   Dimensions do not match.
   */
  class WrongDim : public Error
  {

  public:
    //! Constructor.
    WrongDim (const string& function = "", const string& comment = "") :
        Error (function, comment, "Wrong dimension")
    {
      // Nothing
    }

  };

  ////////////////
  // WRONGINDEX //
  ////////////////

  //! Wrong index.
  /*!
   The index is out of range.
   */
  class WrongIndex : public Error
  {

  public:
    WrongIndex (const string& function = "", const string& comment = "") :
        Error (function, comment, "Index out of range")
    {
      // Nothing
    }

  };

  /////////////
  // IOERROR //
  /////////////

  //! An input/output operation failed.
  class IOError : public Error
  {

  public:
    IOError (const string& function = "", const string& comment = "") :
        Error (function, comment, "An input/output operation failed")
    {
      // Nothing
    }

  };

  ///////////////
  // UNDEFINED //
  ///////////////

  //! Undefined function.
  class Undefined : public Error
  {

  public:
    Undefined (const string& function = "", const string& comment = "") :
        Error ("", comment, Description (function))
    {
      // Nothing
    }

  private:
    //! Returns an error description.
    string
    Description (const string& function)
    {
      if (!function.empty ())
        return string ("Call to undefined function '") + function + "'";
      else
        return "An undefined function was called";
    }

  };

}  // namespace SeldonData.

#define FILE_SELDONDATA_ERRORS_HXX
#endif
