// Copyright (C) 2004-2012, INRIA
// Author(s): Vivien Mallet
//
// This file is part of Talos library, which provides miscellaneous tools to
// make up for C++ lacks and to ease C++ programming.
//
// Talos is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Talos is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Talos. If not, see http://www.gnu.org/licenses/.
//
// For more information, visit the Talos home page:
//     http://vivienmallet.net/lib/talos/


#ifndef TALOS_FILE_STRING_CXX

#include "Date.hxx"
#include "String.hxx"

namespace Talos
{

  //! Default constructor.
  Str::Str()
  {
  }

  //! Copy constructor.
  /*!
    \param[in] s 'Str' instance to be copied.
  */
  Str::Str(const Str& s)
  {
    output_ << s;
  }

  //! Conversion to string.
  Str::operator std::string() const
  {
    return output_.str();
  }

  //! Adds an element to the string.
  /*!
    \param[in] input element added at the end of the string.
  */
  template <class T>
  Str& Str::operator << (const T& input)
  {
    output_ << input;
    return *this;
  }

  //! Adds an element to an instance of 'Str'.
  /*!
    \param[in] s 'Str' instance.
    \param[in] input element added at the end of the string.
  */
  template <class T>
  Str operator + (const Str& s, const T& input)
  {
    string s_input = s;
    Str output;
    output << s_input << input;
    return output;
  }

  //! Converts a 'str' instance to an 'ostream' instance.
  ostream& operator << (ostream& out, Str& in)
  {
    string output = in;
    out << output;
    return out;
  }

  //! Converts a 'str' instance to an 'ostream' instance.
  ostream& operator << (ostream& out, Str in)
  {
    string output = in;
    out << output;
    return out;
  }

  //! Converts most types to string.
  /*!
    \param input variable to be converted.
    \return A string containing 'input'.
  */
  template<typename T>
  string to_str(const T& input)
  {
    ostringstream output;
    output << input;
    return output.str();
  }

  //! Fills a string.
  /*!
    \param input string to be filled.
    \param l (optional) width of the output string. Default: 2.
    \param c (optional) char with which the string will be filled.
    Default: ' '.
    \param flags (optional) format flags. Default: ostringstream::right.
    \return A string containing 'input' (filled).
  */
  string fill(const string& input, int l, char c,
              ostringstream::fmtflags flags)
  {
    std::ostringstream output;
    output.width(l);
    output.fill(c);
    output.flags(flags);
    output << input;
    return output.str();
  }

  //! Converts most types to a filled string.
  /*!
    \param input variable to be converted.
    \param l (optional) width of the output string. Default: 2.
    \param c (optional) char with which the string will be filled.
    Default: ' '.
    \param flags (optional) format flags. Default: ostringstream::right.
    \return A string containing 'input' (filled).
  */
  template<typename T>
  string to_str_fill(const T& input, int l, char c,
                     ostringstream::fmtflags flags)
  {
    ostringstream output;
    output.width(l);
    output.fill(c);
    output.flags(flags);
    output << input;
    return output.str();
  }

  //! Converts string to most types, specially numbers.
  /*!
    \param s string to be converted.
    \param num 's' converted to 'T'.
  */
  template <class T>
  void to_num(const string& s, T& num)
  {
    istringstream str(s);
    str >> num;
  }

  //! Converts string to most types, specially numbers.
  /*!
    \param s string to be converted.
    \return 's' converted to 'T'.
  */
  template <class T>
  T to_num(const string& s)
  {
    T num;
    istringstream str(s);
    str >> num;
    return num;
  }

  //! Converts strings to most types.
  /*!
    \param s string to be converted.
    \param out 's' converted to 'T'.
  */
  template <class T>
  void convert(const string& s, T& out)
  {
    istringstream str(s);
    str >> out;
  }

  //! Sets a string.
  /*!
    \param s input string.
    \param out output string, equal to 's' on exit.
  */
  void convert(const string& s, string& out)
  {
    out = s;
  }

  //! Converts a string to a boolean.
  /*!
    \param s input string.
    \param out output boolean.
  */
  void convert(const string& s, bool& out)
  {
    std::string lower(s);
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   (int(*)(int))tolower);

    if (lower == "true" || lower == "t" || lower == "y" || lower == "yes")
      out = true;
    else if (lower == "false" || lower == "f" || lower == "n"
             || lower == "no")
      out = false;
    else
#ifdef TALOS_DO_NOT_CHECK_BOOLEAN
      {
        istringstream str(s);
        str >> out;
      }
#else
    throw string("Unable to convert \"") + s
      + string("\" to a Boolean. Acceptable strings (case insensitive)")
      + string(" are: true, t, yes, y, false, f, no, n.");
#endif
  }

  //! Converts strings to most types.
  /*!
    \param s input string to be converted.
    \return 's' converted to 'T'.
  */
  template <class T>
  T convert(const string& s)
  {
    T out;
    istringstream str(s);
    str >> out;
    return out;
  }

  //! Converts a string to lower-case string.
  /*!
    \param str string to be converted.
    \return 'str' in lower case.
  */
  string lower_case(string str)
  {
    string lower(str);
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   (int(*)(int))tolower);
    return lower;
  }

  //! Converts a string to upper-case string.
  /*!
    \param str string to be converted.
    \return 'str' in upper case.
  */
  string upper_case(string str)
  {
    string upper(str);
    std::transform(upper.begin(), upper.end(), upper.begin(),
                   (int(*)(int))toupper);
    return upper;
  }

  //! Checks whether a string is a number.
  /*!
    \param str string to be checked.
    \return true if 'str' is a number, false otherwise.
  */
  bool is_num(const string& str)
  {
    if (str == "")
      return false;

    bool mant, mant_a, mant_b, exp;
    string::size_type pos;
    string m, e, m_a, m_b;

    pos = str.find_first_of("eE");
    // Mantissa.
    m = str.substr(0, pos);
    // Exponent.
    e = pos == string::npos ? "" : str.substr(pos + 1);

    exp = pos != string::npos;

    pos = m.find_first_of(".");
    // Mantissa in the form: [m_a].[m_b].
    m_a = m.substr(0, pos);
    // Exponent.
    m_b = pos == string::npos ? "" : m.substr(pos + 1);

    mant = m != "" && m != "-" && m != "+";
    mant_a = m_a != "" && m_a != "-" && m_a != "+";
    mant_b = m_b != "";

    return (mant
            && ((mant_a || mant_b)
                && (!mant_a || is_integer(m_a))
                && (!mant_b || is_unsigned_integer(m_b)))
            && (!exp || is_integer(e)));
  }

  //! Checks whether a string is an integer.
  /*!
    \param str string to be checked.
    \return true if 'str' is an integer, false otherwise.
  */
  bool is_integer(char* str)
  {
    return is_integer(string(str));
  }

  //! Checks whether a string is an integer.
  /*!
    \param str string to be checked.
    \return true if 'str' is an integer, false otherwise.
  */
  bool is_integer(const char* str)
  {
    return is_integer(string(str));
  }

  //! Checks whether a string is an integer.
  /*!
    \param str string to be checked.
    \return true if 'str' is an integer, false otherwise.
  */
  bool is_integer(const string& str)
  {
    bool ans;

    ans = (str.size() > 0 && isdigit(str[0]))
      || (str.size() > 1 && (str[0] == '+' || str[0] == '-'));

    unsigned int i(1);
    while (i < str.size() && ans)
      {
        ans = ans && isdigit(str[i]);
        i++;
      }

    return ans;
  }

  //! Checks whether a string is an unsigned integer.
  /*!
    \param str string to be checked.
    \return true if 'str' is an unsigned integer, false otherwise.
  */
  bool is_unsigned_integer(const string& str)
  {
    bool ans(str.size() > 0);

    unsigned int i(0);
    while (i < str.size() && ans)
      {
        ans = ans && isdigit(str[i]);
        i++;
      }

    return ans;
  }

  //! Checks whether a string is a date.
  /*!
    \param str string to be checked.
    \return true if 'str' is a date, false otherwise.
  */
  bool is_date(const string& str)
  {
    bool ans = true;


    if (!isdigit(str[0]))
      return false;
    try
      {
        Date d(str);
      }
    catch (...)
      {
        ans = false;
      }
    return ans;
  }

  //! Checks whether a string is a valid time interval.
  /*!
    \param str string to be checked.
    \return true if 'str' is a valid time interval, false otherwise.
  */
  bool is_delta(const string& str)
  {
    if ((str.find('h', 0) == string::npos &&
         str.find('d', 0) == string::npos) ||
        (str.find('d', 0) != string::npos &&
         str.find('d', 0) >= str.find('h', 0)))
      return false;

    vector<string> period = split(str, "dh-_");

    if (period.size() > 2 || period.size() == 0)
      return false;

    if (period.size() == 2)
      return (is_integer(period[0]) && is_integer(period[1]));
    else
      return is_integer(period[0]);
  }

  //! Finds and replace a substring.
  /*!
    \param str base string.
    \param old_str substring to be replaced.
    \param new_str substring to be put in place of 'old_str'.
    \return 'str' where 'old_str' was replaced by 'new'str'.
  */
  string find_replace(string str, string old_str, string new_str)
  {
    string::size_type index = str.find(old_str);

    while (index != string::npos)
      {
        str.replace(index, old_str.size(), new_str);
        index = str.find(old_str, index + new_str.size());
      }

    return str;
  }

  //! Trims off a string.
  /*!
    Removes delimiters at each edge of the string.
    \param str string to be trimmed off.
    \param delimiters characters to be removed.
    \return 'str' trimmed off.
  */
  string trim(string str, string delimiters)
  {
    string::size_type index_end = str.find_last_not_of(delimiters);
    string::size_type index_beg = str.find_first_not_of(delimiters);

    if (index_beg == string::npos)
      return "";

    return str.substr(index_beg, index_end - index_beg + 1);
  }

  //! Trims off a string.
  /*!
    Removes delimiters at the beginning of the string.
    \param str string to be trimmed off.
    \param delimiters characters to be removed.
    \return 'str' trimmed off at the beginning.
  */
  string trim_beg(string str, string delimiters)
  {
    string::size_type index = str.find_first_not_of(delimiters);

    if (index == string::npos)
      return "";

    return str.substr(index);
  }

  //! Trims off a string.
  /*!
    Removes delimiters at the end of the string.
    \param str string to be trimmed off.
    \param delimiters characters to be removed.
    \return 'str' trimmed off at the end.
  */
  string trim_end(string str, string delimiters)
  {
    string::size_type index = str.find_last_not_of(delimiters);

    if (index == string::npos)
      return "";

    return str.substr(0, index + 1);
  }

  //! Splits a string.
  /*!
    The string is split according to delimiters and elements are stored
    in the vector 'vect'.
    \param str string to be split.
    \param vect (output) vector containing elements of the string.
    \param delimiters (optional) delimiters. Default: " \n\t".
  */
  template <class T>
  void split(string str, vector<T>& vect, string delimiters)
  {
    vect.clear();

    T tmp;
    string::size_type index_beg, index_end;

    index_beg = str.find_first_not_of(delimiters);

    while (index_beg != string::npos)
      {
        index_end = str.find_first_of(delimiters, index_beg);
        convert(str.substr(index_beg, index_end == string::npos ?
                           string::npos : (index_end - index_beg)), tmp);
        vect.push_back(tmp);
        index_beg = str.find_first_not_of(delimiters, index_end);
      }
  }

  //! Splits a string.
  /*!
    The string is split according to delimiters.
    \param str string to be split.
    \param delimiters (optional) delimiters. Default: " \n\t".
    \return A vector containing elements of the string.
  */
  vector<string> split(string str, string delimiters)
  {
    vector<string> vect;
    split(str, vect, delimiters);
    return vect;
  }

  //! Extracts markups from a string.
  /*!
    The string is split into markups and elements.
    \param str string to be split.
    \param elements (output) vector containing markups (without their tags)
    and elements of the string.
    \param is_markup (output) booleans set to true for each markup found in
    'elements'.
    \param delimiters (optional) markup tags. Default: "$".
    \note A markup is a field delimited by two tags.
  */
  template <class T>
  void split_markup(string str, vector<T>& elements,
                    vector<bool>& is_markup, string delimiters)
  {
    elements.clear();
    is_markup.clear();

    string current;
    T tmp;
    bool is_mark;
    string::size_type index_beg, index_end, markup_end;

    index_beg = 0;

    while (index_beg != string::npos && index_beg < str.size())
      {
        is_mark = false;
        // Stops at the first delimiter.
        index_end = str.find_first_of(delimiters, index_beg);

        // No more delimiter.
        if (index_end == string::npos)
          {
            current = str.substr(index_beg, string::npos);
            markup_end = string::npos;
          }
        // One character left, added later.
        else if (index_end == str.size() - 1)
          {
            current = str.substr(index_beg, index_end - index_beg);
            markup_end = string::npos;
          }
        // There may have a markup there. Searches for the end of this markup.
        else
          {
            current = str.substr(index_beg, index_end - index_beg);
            markup_end = str.find_first_of(delimiters, index_end + 1);
          }

        // No markup: e.g. "$$" means "$".
        if (markup_end != string::npos && markup_end - 1 == index_end)
          {
            current += str[markup_end];
            index_end = markup_end;
          }
        // A markup is detected.
        else if (markup_end != string::npos)
          {
            // The previous element, which is not a markup, is added.
            if (current != "")
              if (elements.size() > 0 && !is_markup[elements.size() - 1])
                elements[elements.size() - 1] += current;
              else
                {
                  elements.push_back(current);
                  is_markup.push_back(false);
                }
            // The markup is extracted.
            current = str.substr(index_end, markup_end - index_end + 1);
            index_end = markup_end;
            is_mark = true;
          }
        // Gets the rest of the string.
        else if (index_end != string::npos)
          current += str.substr(index_end, markup_end);

        // If not a markup, then adds the string to the last stored element
        // (if it is not a markup itself).
        if (elements.size() > 0 && !is_markup[elements.size() - 1]
            && !is_mark)
          elements[elements.size() - 1] += current;
        else
          {
            // Removes markup tags, for markups.
            if (is_mark)
              current = trim(current, delimiters);
            elements.push_back(current);
            is_markup.push_back(is_mark);
          }

        // Moves forward.
        index_beg = markup_end == string::npos ? string::npos : index_end + 1;

      }
  }

  //! Displays a vector.
  /*!
    \param v the vector to de displayed.
  */
  template <class T>
  void print(const vector<T>& v)
  {
    if (v.size() != 0)
      cout << v[0];
    for (unsigned int i = 1; i < v.size(); i++)
      cout << '\t' << v[i];
    cout << endl;
  }

}  // namespace Talos.


#define TALOS_FILE_STRING_CXX
#endif
