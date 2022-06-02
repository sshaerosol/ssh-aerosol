// Copyright (C) 2004-2007, INRIA
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


#ifndef TALOS_FILE_FILES_CXX

#include "String.hxx"
#include "Files.hxx"

#include <algorithm>
#include <map>
#include <set>

namespace Talos
{

  //! Tests whether a file exists.
  /*!
    \param file_name file name.
    \return true if the file exists, false otherwise.
    \note This function returns false if the user has not the rights to
    read the file, even if the file exists.
  */
  bool exists(string file_name)
  {
    ifstream file_stream(file_name.c_str(), ifstream::in);
    bool ans = file_stream.is_open();
    file_stream.close();

    return ans;
  }

  //! Returns a file size.
  /*!
    \param file_name file name.
    \return The file size in bytes.
  */
  unsigned long file_size(string file_name)
  {
    ifstream file_stream(file_name.c_str(), ifstream::in);
    streampos position = file_stream.tellg();
    file_stream.seekg(0, ios::end);

    return (file_stream.tellg() - position);
  }

  //! Returns a stream size.
  /*!
    \param stream the stream.
    \return The stream size in bytes.
  */
  unsigned long stream_size(istream& stream)
  {
    streampos position = stream.tellg();
    stream.seekg(0, ios::end);
    unsigned long length = stream.tellg() - position;
    stream.seekg(position, ios::beg);

    return length;
  }

  //! Checks whether a stream is empty.
  /*!
    \param stream the stream.
    \return true if the stream is empty, false otherwise.
  */
  bool is_emptystream(istream& stream)
  {
    streampos position = stream.tellg();
    istream::iostate state = stream.rdstate();

    stream.get();
    bool res = !stream.good();

    stream.clear(state);
    stream.seekg(position);

    return res;
  }

  //! Checks whether a stream contains an element.
  /*!
    Checks whether a stream contains an element that may be extracted
    through 'operator >>'.
    \param stream the stream.
    \return true if the stream has an element, false otherwise.
  */
  bool has_element(istream& stream)
  {
    streampos position = stream.tellg();
    istream::iostate state = stream.rdstate();

    string element;
    bool res = bool(stream >> element); // YK

    stream.clear(state);
    stream.seekg(position);

    return res;
  }

  //! Checks whether a numerical value satisfies a list of constraints.
  /*!
    \param value the numerical value.
    \param constraint the list of constraints. The constraints are delimited
    by |. The supported constraints are: positive, strictly positive,
    negative, strictly negative, non zero, integer, > x, >= x, < x, <= x, != x
    y z, = x y z.
    \return true if the constraints are satisfied, false otherwise.
  */
  template <class T>
  bool satisfies_constraint(T value, string constraint)
  {
    vector<string> constraint_list = split(constraint, "|");

    string expression, str;
    T val;
    for (int i = 0; i < int(constraint_list.size()); i++)
      {
        expression = trim(constraint_list[i]);
        if (expression.size() < 2)
          throw "Error in satisfies_constraint: the constraint \""
            + expression + "\" is not supported.";
        if (expression[0] == '<')
          if (expression[1] == '=')
            {
              str = trim(expression.substr(2));
              if (!is_num(str))
                throw "Error in satisfies_constraint: the constraint \""
                  + expression + "\" cannot be parsed or is not supported.";
              to_num(str, val);
              if (value > val)
                return false;
            }
          else
            {
              str = trim(expression.substr(1));
              if (!is_num(str))
                throw "Error in satisfies_constraint: the constraint \""
                  + expression + "\" cannot be parsed or is not supported.";
              to_num(str, val);
              if (value >= val)
                return false;
            }
        else if (expression[0] == '>')
          if (expression[1] == '=')
            {
              str = trim(expression.substr(2));
              if (!is_num(str))
                throw "Error in satisfies_constraint: the constraint \""
                  + expression + "\" cannot be parsed or is not supported.";
              to_num(str, val);
              if (value < val)
                return false;
            }
          else
            {
              str = trim(expression.substr(1));
              if (!is_num(str))
                throw "Error in satisfies_constraint: the constraint \""
                  + expression + "\" cannot be parsed or is not supported.";
              to_num(str, val);
              if (value <= val)
                return false;
            }
        else if (expression.substr(0, 2) == "!=")
          {
            str = trim(expression.substr(2));
            vector<string> number = split(str);
            for (int j = 0; j < int(number.size()); j++)
              {
                if (!is_num(number[j]))
                  throw "Error in satisfies_constraint: the constraint \""
                    + expression + "\" cannot be parsed or is not supported.";
                to_num(number[j], val);
                if (value == val)
                  return false;
              }
          }
        else if (expression[0] == '=')
          {
            str = trim(expression.substr(1));
            vector<string> number = split(str);
            bool acceptable = false;
            for (int j = 0; j < int(number.size()); j++)
              {
                if (!is_num(number[j]))
                  throw "Error in satisfies_constraint: the constraint \""
                    + expression + "\" cannot be parsed or is not supported.";
                to_num(number[j], val);
                acceptable = acceptable || value == val;
              }
            if (!acceptable)
              return false;
          }
        else if (expression == "positive")
          {
            if (value < T(0))
              return false;
          }
        else if (expression == "strictly positive")
          {
            if (value <= T(0))
              return false;
          }
        else if (expression == "negative")
          {
            if (value > T(0))
              return false;
          }
        else if (expression == "strictly negative")
          {
            if (value >= T(0))
              return false;
          }
        else if (expression == "non zero")
          {
            if (value == T(0))
              return false;
          }
        else if (expression == "integer")
          {
            if (value != T(int(value)))
              return false;
          }
        else
          throw "Error in satisfies_constraint: the constraint \""
            + expression + "\" cannot be parsed or is not supported.";
      }

    return true;
  }

  //! Formats a list of constraints in a string readable for human beings.
  /*!
    \param constraint the list of constraints. The constraints are delimited
    by |. The supported constraints are: positive, strictly positive,
    negative, strictly negative, non zero, integer, > x, >= x, < x, <= x, != x
    y z, = x y z.
    \return The constraints in a readable string.
  */
  string show_constraint(string constraint)
  {
    vector<string> constraint_list = split(constraint, "|");

    string output = "";

    string expression, str, termination;
    for (int i = 0; i < int(constraint_list.size()); i++)
      {
        if (i != int(constraint_list.size()) - 1)
          termination = ";\n";
        else
          termination = ".";
        expression = trim(constraint_list[i]);
        if (expression.size() < 2)
          throw "Error in show_constraint: the constraint \""
            + expression + "\" cannot be parsed.";
        if (expression[0] == '<')
          if (expression[1] == '=')
            {
              str = trim(expression.substr(2));
              if (!is_num(str))
                throw "Error in show_constraint: the constraint \""
                  + expression + "\" cannot be parsed or is not supported.";
              output += " - Value less than " + str + termination;
            }
          else
            {
              str = trim(expression.substr(1));
              if (!is_num(str))
                throw "Error in show_constraint: the constraint \""
                  + expression + "\" cannot be parsed or is not supported.";
              output += " - Value strictly less than " + str + termination;
            }
        else if (expression[0] == '>')
          if (expression[1] == '=')
            {
              str = trim(expression.substr(2));
              if (!is_num(str))
                throw "Error in show_constraint: the constraint \""
                  + expression + "\" cannot be parsed or is not supported.";
              output += " - Value greater than " + str + termination;
            }
          else
            {
              str = trim(expression.substr(1));
              if (!is_num(str))
                throw "Error in show_constraint: the constraint \""
                  + expression + "\" cannot be parsed or is not supported.";
              output += " - Value strictly greater than " + str + termination;
            }
        else if (expression.substr(0, 2) == "!=")
          {
            str = trim(expression.substr(2));
            vector<string> number = split(str);
            output += " - Value different from ";
            for (int j = 0; j < int(number.size()); j++)
              {
                if (!is_num(number[j]))
                  throw "Error in show_constraint: the constraint \""
                    + expression + "\" cannot be parsed or is not supported.";
                output += number[j];
                if (j == int(number.size()) - 2)
                  output += " and ";
                else if (j != int(number.size()) - 1)
                  output += ", ";
              }
            output += termination;
          }
        else if (expression[0] == '=')
          {
            str = trim(expression.substr(1));
            vector<string> number = split(str);
            output += " - Value equal to ";
            for (int j = 0; j < int(number.size()); j++)
              {
                if (!is_num(number[j]))
                  throw "Error in show_constraint: the constraint \""
                    + expression + "\" cannot be parsed or is not supported.";
                output += number[j];
                if (j == int(number.size()) - 2)
                  output += " or ";
                else if (j != int(number.size()) - 1)
                  output += ", ";
              }
            output += termination;
          }
        else if (expression == "positive")
          output += " - Positive value" + termination;
        else if (expression == "strictly positive")
          output += " - Strictly positive value" + termination;
        else if (expression == "negative")
          output += " - Negative value" + termination;
        else if (expression == "strictly negative")
          output += " - Strictly negative value" + termination;
        else if (expression == "non zero")
          output += " - Non-zero value" + termination;
        else if (expression == "integer")
          output += " - Integral value" + termination;
        else
          throw "Error in show_constraint: the constraint \""
            + expression + "\" cannot be parsed or is not supported.";
      }

    return output;
  }


  /////////////////
  // SEARCHSCOPE //
  /////////////////

  //! A scope opened when searching for a field.
  class SearchScope
  {
  public:
    //! Main constructor.
    /*! Begins the search scope.
      \param stream the stream where the search takes place.
      \param searching the term being searched for.
    */
    SearchScope(ExtStream& stream, const string& searching)
      : searching_(stream.searching_)
    {
      searching_ = searching;
#ifdef TALOS_DEBUG
      AddToRegister(stream);
#endif
    }

    //! Main constructor.
    /*! Begins the search scope for multiple streams.
     */
    SearchScope(ConfigStreams& multi_stream, const string& searching)
      : searching_(multi_stream.searching_)
    {
      searching_ = searching;

      // No real scope for sub-streams, just registering.
      for (int i = 0; i < (int) multi_stream.streams_.size(); ++i)
        SearchScope(*multi_stream.streams_[i], searching);
    }

    //! Destructor.
    /*! Ends the search scope.
     */
    ~SearchScope()
    {
      searching_ = "";
    }

  private:
    SearchScope(const SearchScope&);
    string& searching_;

#ifdef TALOS_DEBUG
    void AddToRegister(ExtStream& stream)
    {
      reg.file[stream.file_name_].insert(searching_);
      string& delimiter = reg.delimiter[stream.file_name_];
      if (delimiter.empty())
        delimiter = stream.delimiters_;
      else if (delimiter != stream.delimiters_)
        {
          throw "\"" + stream.file_name_ + "\" has been opened with different "
            "delimiters: once with \"" + delimiter + "\", "
            "another with \"" + stream.delimiters_ + "\"";
        }
    }

    //! Fields that has been searched for through the configuration files.
    static struct Register
    {
      typedef std::set<string> FieldSet;
      typedef map<string, FieldSet> FileRegister;
      typedef map<string, string> DelimiterRegister;
      FileRegister file;
      DelimiterRegister delimiter;

      //! Destructor.
      //! Checks that all configuration field names has been requested.
      ~Register()
      {
        bool has_warning = false;
        for (FileRegister::iterator it = file.begin(); it != file.end(); ++it)
          {
            const string& filename = it->first;
            const FieldSet& searched_field = it->second;
            const string& delimiters = delimiter[filename];

            ConfigStream cfg(filename);
            FieldSet field_list;

            // Adds the sections to the field list.
            string element;
            while (cfg.GetRawElement(element))
              if (cfg.IsSection(element))
                field_list.insert(element);

            // Adds variables to the field list.
            cfg.Rewind();
            string line;
            vector<string> member_list;
            vector<string> word_list;
            while (cfg.ExtStream::GetLine(line))
              {
                vector<string> variable_list;

                split(line, member_list, ":=");
                int member_count = (int) member_list.size();
                for (int i = 0; i < member_count; i += 2)
                  {
                    const string& token = member_list[i];

                    split(token, word_list, delimiters);
                    int word_count = word_list.size();
                    if (member_count == 1  // no affectation symbol
                        && word_count > 2) // and not a pair of words
                      break;               // => a line of raw values

                    string variable_name = word_list.back();
                    if (is_num(variable_name))
                      {
                        variable_list.clear();
                        break; // This was actually a line of raw values.
                      }
                    variable_list.push_back(variable_name);
                  }

                for (int i = 0; i < (int) variable_list.size(); ++i)
                  field_list.insert(variable_list[i]);
              }

            FieldSet unused_field;
            set_difference(field_list.begin(), field_list.end(),
                           searched_field.begin(), searched_field.end(),
                           inserter(unused_field, unused_field.begin()));

            // Checks that every fields were searched for.
            if (!unused_field.empty())
              {
                has_warning = true;
                cerr << "[INFO] === in \"" << filename << "\" ===" << endl;
              }
            for (FieldSet::iterator it = unused_field.begin();
                 it != unused_field.end(); ++it)
              cerr << "[INFO]   '"
                   << *it << "' field was never used " << endl;
          }
        if (has_warning)
          cerr << "[INFO] Caveat of the unused field detection:\n"
            " - A field name with multiple occurrences can be unused for some\n"
            "   of them without being listed here.\n"
            " - Some unused fields might actually be used by another reader\n"
            "   than the Talos config reader." << endl;
      }
    } reg;
#endif
  };

#ifdef TALOS_DEBUG
  SearchScope::Register SearchScope::reg = SearchScope::Register();
#endif

  ///////////////
  // EXTSTREAM //
  ///////////////

  //! Default constructor.
  /*! Nothing is performed.
   */
  ExtStream::ExtStream():
    comments_("#%"), delimiters_(" \t:=|\n,;\r\x0D\x0A"), searching_("")
  {
  }

  //! Main constructor.
  /*! Opens a file.
    \param file_name file to be opened.
  */
  ExtStream::ExtStream(string file_name,
                       string comments,
                       string delimiters):
    ifstream(file_name.c_str(), ifstream::binary), file_name_(file_name),
    comments_(comments), delimiters_(delimiters), searching_("")
  {
    if (!this->is_open())
      throw string("Unable to open file \"") + file_name + "\".";
  }

  //! Destructor.
  /*! Closes the stream.
   */
  ExtStream::~ExtStream()
  {
    this->close();
  }

  //! Checks whether a line should be discarded.
  /*!
    \param line line to be checked.
    \return true if the line should be discarded, false otherwise.
  */
  bool ExtStream::Discard(string line) const
  {
    size_t first = line.find_first_not_of(delimiters_);
    return ((first == string::npos)
            || (comments_.find_first_of(line[first]) != string::npos));
  }

  //! Skips discarded lines.
  /*!
    Extracts discarded lines.
    \return A reference to the current stream.
  */
  ExtStream& ExtStream::SkipDiscarded()
  {
    std::streampos position;
    while ((!is_emptystream(*this)) && (Discard(PeekFullLine(position))))
      this->seekg(position);
    return *this;
  }

  //! Sets the characters considered as delimiters.
  /*!
    \param delimiters delimiters.
  */
  void ExtStream::SetDelimiters(string delimiters)
  {
    delimiters_ = delimiters;
  }

  //! Sets the characters that denote a comment line.
  /*!
    \param comments the characters that denote a comment line.
  */
  void ExtStream::SetComments(string comments)
  {
    comments_ = comments;
  }

  //! Returns the characters considered as delimiters..
  /*!
    \return Delimiters.
  */
  string ExtStream::GetDelimiters() const
  {
    return delimiters_;
  }

  //! Returns the characters that denote a comment line.
  /*!
    \return The characters that denote a comment line.
  */
  string ExtStream::GetComments() const
  {
    return comments_;
  }

  //! Returns the name of the file that was opened.
  /*!
    \return The name of the file that was opened.
  */
  string ExtStream::GetFileName() const
  {
    return file_name_;
  }

  //! Skips delimiters.
  /*!
    Extracts following delimiters from the string, until another character
    is found.
    \return A reference to the current stream.
  */
  ExtStream& ExtStream::SkipDelimiters()
  {
    while (this->good()
           && (delimiters_.find_first_of(char(this->peek()))
               != string::npos))
      this->get();
    return *this;
  }

  //! Removes delimiters at both ends of a string.
  /*!
    Removes delimiters at the beginning and at the end of a string.
    \param str string.
    \return The string without delimiters at both ends.
  */
  string ExtStream::RemoveDelimiters(const string& str) const
  {
    string::size_type index = str.find_first_not_of(delimiters_);
    if (index == string::npos)
      return "";
    return str.substr(index,
                      str.find_last_not_of(delimiters_) - index + 1);
  }

  //! Skips discarded lines and delimiters.
  /*!
    Extracts discarded lines and delimiters.
    \return A reference to the current stream.
  */
  ExtStream& ExtStream::Skip()
  {
    this->SkipDiscarded();
    return this->SkipDelimiters();
  }

  //! Opens a file.
  /*!
    \param file_name file name.
    \param mode (optional) flags describing the requested I/O mode for the
    file.  Default: in.
    \note If a file was previously opened, it is closed and the stream is
    cleared.
  */
  void ExtStream::Open(string file_name, openmode mode)
  {
    this->close();
    this->clear();
    this->open(file_name.c_str(), mode);

    file_name_ = file_name;

    if (!this->is_open())
      throw string("Unable to open file \"") + file_name + "\".";
  }

  //! Closes the current file.
  /*!
    \note The stream is cleared.
  */
  void ExtStream::Close()
  {
    this->close();
    this->clear();

    file_name_ = "";
  }

  //! Checks whether the stream is empty.
  /*!
    Checks whether the stream has still valid elements to be read.
    \return 'true' is the stream is empty, 'false' otherwise.
  */
  bool ExtStream::IsEmpty()
  {
    string tmp;
    return !this->PeekElement(tmp);
  }

  //! Rewinds the stream.
  /*!
    Goes back to the beginning of the stream and clears the control state.
    \return A reference to the current stream.
  */
  ExtStream& ExtStream::Rewind()
  {
    this->clear();
    this->seekg(0, ifstream::beg);

    return *this;
  }

  //! Returns the next line.
  /*!
    \return The next line.
  */
  string ExtStream::GetFullLine()
  {
    string line;

    std::getline(*this, line);

    return line;
  }

  //! Returns the next line.
  /*!
    \param line (output) the next line.
  */
  bool ExtStream::GetFullLine(string& line)
  {
    return bool(std::getline(*this, line)); // YK
  }

  //! Returns the next line without extracting it from the stream.
  /*!
    \return The next line.
  */
  string ExtStream::PeekFullLine()
  {
    std::streampos position = this->tellg(); 
    iostate state = this->rdstate();

    string line;
    std::getline(*this, line);

    this->clear(state);
    this->seekg(position);

    return line;
  }

  //! Returns the next line without extracting it from the stream.
  /*!
    \param position (output) the position of the line following the next line.
    \return The next line.
  */
  string ExtStream::PeekFullLine(std::streampos& position)
  {
    std::streampos position_back = this->tellg();
    iostate state = this->rdstate();

    string line;
    std::getline(*this, line);

    position = this->tellg();

    this->clear(state);
    this->seekg(position_back);

    return line;
  }

  //! Returns the next line without extracting it from the stream.
  /*!
    \param line (output) the next line.
    \return true if a line has been found, false otherwise.
  */
  bool ExtStream::PeekFullLine(string& line)
  {
    std::streampos position = this->tellg();
    iostate state = this->rdstate();

    bool success = bool(std::getline(*this, line));

    this->clear(state);
    this->seekg(position);

    return success;
  }

  //! Skips full lines.
  /*!
    \param nb number of lines to be skipped.
  */
  void ExtStream::SkipFullLines(int nb)
  {
    for (int i = 0; i < nb; i++)
      this->GetFullLine();
  }

  //! Returns the next valid line.
  /*!
    Returns the next valid line, i.e. the next line that is
    not a line to be discarded and from which comments have been extracted.
    \return The next valid line.
  */
  string ExtStream::GetLine()
  {
    bool not_end;
    string line;
    string::size_type index(0), index_tmp;

    this->Skip();
    line = GetFullLine();

    while ((not_end = ((index_tmp
                        = line.substr(index).find_first_of(comments_))
                       != string::npos))
           && (delimiters_.find_first_of(line[(index += index_tmp) - 1])
               == string::npos)
           && (not_end = (++index != line.size())));

    if (not_end)
      index --;
    else
      index = line.size();

    while (delimiters_.find_first_of(line[--index]) != string::npos);

    return line.substr(0, index + 1);
  }

  //! Returns the next valid line.
  /*!
    Returns the next valid line, i.e. the next line that is
    not a line to be discarded and from which comments have been extracted.
    \param line (output) the next valid line.
  */
  bool ExtStream::GetLine(string& line)
  {
    bool not_end, success;
    string::size_type index(0), index_tmp;

    this->Skip();
    success = GetFullLine(line);

    while ((not_end = ((index_tmp
                        = line.substr(index).find_first_of(comments_))
                       != string::npos))
           && (delimiters_.find_first_of(line[(index += index_tmp) - 1])
               == string::npos)
           && (not_end = (++index != line.size())));

    if (not_end)
      index --;
    else
      index = line.size();

    while (delimiters_.find_first_of(line[--index]) != string::npos);

    line = line.substr(0, index + 1);

    return success;
  }

  //! Returns the next valid line without extracting it from the stream.
  /*!
    Returns the next valid line, i.e. the next line that is
    not a line to be discarded and from which comments have been extracted.
    Nothing is extracted from the stream.
    \return The next valid line.
  */
  string ExtStream::PeekLine()
  {
    string line;

    std::streampos initial_position = this->tellg();
    iostate state = this->rdstate();

    line = this->GetLine();

    this->clear(state);
    this->seekg(initial_position);

    return line;
  }

  //! Returns the next valid line without extracting it from the stream.
  /*!
    Returns the next valid line, i.e. the next line that is
    not a line to be discarded and from which comments have been extracted.
    Nothing is extracted from the stream.
    \param position (output) the position of the line following the next valid
    line.
    \return The valid line.
  */
  string ExtStream::PeekLine(std::streampos& position)
  {
    std::streampos position_back = this->tellg();
    iostate state = this->rdstate();

    string line = this->GetLine();

    position = this->tellg();

    this->clear(state);
    this->seekg(position_back);

    return line;
  }

  //! Returns the next valid line without extracting it from the stream.
  /*!
    Returns the next valid line, i.e. the next line that is
    not a line to be discarded and from which comments have been extracted.
    Nothing is extracted from the stream.
    \param line (output) the next valid line.
  */
  bool ExtStream::PeekLine(string& line)
  {
    std::streampos initial_position = this->tellg();
    iostate state = this->rdstate();

    bool success = this->GetLine(line);

    this->clear(state);
    this->seekg(initial_position);

    return success;
  }

  //! Skips valid lines.
  /*!
    \param nb number of lines to be skipped.
  */
  void ExtStream::SkipLines(int nb)
  {
    for (int i = 0; i < nb; i++)
      this->GetLine();
  }

  //! Sets the position of the get pointer after a given element.
  /*!
    Sets the position of the get pointer exactly after a given element.
    \param element the element to be found.
    \return true if the element was found, false otherwise.
  */
  bool ExtStream::Find(string element)
  {
    SearchScope s(*this, element);

    string elt;
    while (GetElement(elt) && elt != element);

    if (elt == "")
      throw string("Error in ExtStream::Find: \"")
        + element + string("\" not found in \"") + file_name_ + "\".";

    return elt == element;
  }

  //! Sets the position of the get pointer after a given element.
  /*!
    Sets the position of the get pointer exactly after a given element.
    \param element the element to be found from the beginning of the stream.
    \return true if the element was found, false otherwise.
  */
  bool ExtStream::FindFromBeginning(string element)
  {
    this->Rewind();
    return this->Find(element);
  }

  //! Returns the next valid element.
  /*!
    Returns the next valid element, i.e. the next element that is
    not in a line to be discarded.
    \return The next valid element.
  */
  string ExtStream::GetElement()
  {
    std::streampos position;
    string element;
    string::size_type index, length;

    while ((!is_emptystream(*this)) && (Discard(PeekFullLine(position))))
      this->seekg(position);
    element = PeekFullLine();

    index = element.find_first_not_of(delimiters_);
    if (index != string::npos)
      {
        length = element.find_first_of(delimiters_, index);
        length = length == string::npos ? element.size() - index
          : length - index;
        element = element.substr(index, length);
      }
    else
      {
        index = element.size();
        length = 0;
        element = "";
      }

    this->seekg(index + length, ifstream::cur);

    return element;
  }

  //! Gets the next valid element.
  /*!
    Gets the next valid element, i.e. the next element that is
    not in a line to be discarded.
    \param element (output) the next valid element.
  */
  template <class T>
  bool ExtStream::GetElement(T& element)
  {
    string str = GetElement();
    convert(str, element);

    return (str != "");
  }

  //! Gets the next valid element.
  /*!
    Gets the next valid element, i.e. the next element that is
    not in a line to be discarded.
    \param element (output) the next valid element.
  */
  template <class T>
  bool ExtStream::GetRawElement(T& element)
  {
    string str = ExtStream::GetElement();
    convert(str, element);

    return (str != "");
  }

  //! Returns the next valid element without extracting it from the stream.
  /*!
    Returns the next valid element, i.e. the next element that is
    not in a line to be discarded.
    \return The next valid element.
  */
  string ExtStream::PeekElement()
  {
    std::streampos initial_position;
    string element;

    initial_position = this->tellg();
    iostate state = this->rdstate();

    element = GetElement();

    this->clear(state);
    this->seekg(initial_position);

    return element;
  }

  //! Gets the next valid element without extracting it from the stream.
  /*!
    Gets the next valid element, i.e. the next element that is
    not in a line to be discarded.
    \param element (output) the next valid element.
  */
  template <class T>
  bool ExtStream::PeekElement(T& element)
  {
    std::streampos initial_position;
    bool success;

    initial_position = this->tellg();
    iostate state = this->rdstate();

    success = GetElement(element);

    this->clear(state);
    this->seekg(initial_position);

    return success;
  }

  //! Skips valid elements.
  /*!
    \param nb number of valid elements to be skipped.
  */
  void ExtStream::SkipElements(int nb)
  {
    for (int i = 0; i < nb; i++)
      GetElement();
  }

  //! Returns the next valid number.
  /*!
    Returns the next valid number, i.e. the next number that is
    not in a line to be discarded.
    \return The next valid number.
  */
  double ExtStream::GetNumber()
  {
    string element;
    while (GetElement(element) && !is_num(element));

    return is_num(element) ? to_num<double>(element) : 0.;
  }

  //! Gets the next valid number.
  /*!
    Gets the next valid number, i.e. the next number that is
    not in a line to be discarded.
    \param element (output) the next valid number.
  */
  template <class T>
  bool ExtStream::GetNumber(T& number)
  {
    string element;
    bool success;
    while ((success = GetElement(element)) && !is_num(element));

    number = is_num(element) ? to_num<T>(element) : T(0);

    return success;
  }

  //! Returns the next valid number without extracting it from the stream.
  /*!
    Returns the next valid number, i.e. the next number that is
    not in a line to be discarded.
    \return The next valid number.
  */
  double ExtStream::PeekNumber()
  {
    std::streampos initial_position = this->tellg();
    iostate state = this->rdstate();

    string element;
    while (GetElement(element) && !is_num(element));

    this->clear(state);
    this->seekg(initial_position);

    return is_num(element) ? to_num<double>(element) : 0.;
  }

  //! Gets the next valid number without extracting it from the stream.
  /*!
    Gets the next valid number, i.e. the next number that is
    not in a line to be discarded.
    \param number (output) the next valid number.
  */
  template <class T>
  bool ExtStream::PeekNumber(T& number)
  {
    std::streampos initial_position = this->tellg();
    iostate state = this->rdstate();

    string element;
    bool success;
    while ((success = GetElement(element)) && !is_num(element));

    number = is_num(element) ? to_num<T>(element) : T(0);

    this->clear(state);
    this->seekg(initial_position);

    return success;
  }

  //! Skips numbers.
  /*!
    \param nb number of numbers to be skipped.
  */
  void ExtStream::SkipNumbers(int nb)
  {
    for (int i = 0; i < nb; i++)
      this->GetNumber();
  }

  //! Gets the value of a given variable.
  /*!
    Gets the value of a given variable, i.e. the next valid
    (not in a discarded line) element following the variable name.
    \param name the name of the variable.
    \return the value of the variable.
  */
  string ExtStream::GetValue(string name)
  {
    SearchScope s(*this, name);

    string element;
    while (GetElement(element) && element != name);

    if (element != name)
      throw string("Error in ExtStream::GetValue: \"")
        + name + string("\" not found in \"") + file_name_ + "\".";

    return GetElement();
  }

  //! Gets the value of a given variable without extracting from the stream.
  /*!
    Gets the value of a given variable, i.e. the next valid
    (not in a discarded line) element following the variable name.
    Nothing is extracted from the stream.
    \param name the name of the variable.
    \return the value associated with the variable.
  */
  string ExtStream::PeekValue(string name)
  {
    std::streampos initial_position = this->tellg();
    iostate state = this->rdstate();

    string element = this->GetValue(name);

    this->clear(state);
    this->seekg(initial_position);

    return element;
  }

  //! Gets the value of a given variable.
  /*!
    Gets the (numerical) value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    \param name the name of the variable.
    \param value value associated with the variable.
  */
  template <class T>
  void ExtStream::GetValue(string name, T& value)
  {
    SearchScope s(*this, name);

    string element;
    while (GetElement(element) && element != name);

    if (element != name)
      throw string("Error in ExtStream::GetValue: \"")
        + name + string("\" not found in \"") + file_name_ + "\".";
    if (!this->GetElement(element))
      throw string("Error in ExtStream::GetValue: unable to read value of \"")
        + name + string("\" in \"") + file_name_ + "\".";
    if (!is_num(element))
      throw string("Error in ExtStream::GetValue: the value of \"") + name
        + string("\" in \"") + file_name_ + string("\" is \"") + element
        + "\", but it should be a number.";

    value = to_num<T>(element);
  }

  //! Gets the value of a given variable.
  /*!
    Gets the integral value of a given variable, i.e. the next valid (not in a
    discarded line) number or element following the variable name.
    \param name the name of the variable.
    \param value value associated with the variable.
  */
  void ExtStream::GetValue(string name, int& value)
  {
    SearchScope s(*this, name);

    string element;
    while (GetElement(element) && element != name);

    if (element != name)
      throw string("Error in ExtStream::GetValue: \"")
        + name + string("\" not found in \"") + file_name_ + "\".";
    if (!this->GetElement(element))
      throw string("Error in ExtStream::GetValue: unable to read value of \"")
        + name + string("\" in \"") + file_name_ + "\".";
    if (!is_integer(element))
      throw string("Error in ExtStream::GetValue: the value of \"") + name
        + string("\" in \"") + file_name_ + string("\" is \"") + element
        + "\", but it should be an integer.";

    value = to_num<int>(element);
  }

  /*! \brief Gets the value of a given variable without extracting them from
    the stream. */
  /*!
    Gets the (numerical) value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    Nothing is extracted from the stream.
    \param name the name of the variable.
    \param value value associated with the variable.
  */
  template <class T>
  void ExtStream::PeekValue(string name, T& value)
  {
    std::streampos initial_position = this->tellg();
    iostate state = this->rdstate();

    this->GetValue(name, value);

    this->clear(state);
    this->seekg(initial_position);
  }

  //! Gets the value of a given variable.
  /*!
    Gets the (numerical) value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    \param name the name of the variable.
    \param min the minimum value that the variable should take.
    \param max the maximum value that the variable should take.
    \param value value associated with the variable.
  */
  template <class T>
  void ExtStream::GetValue(string name, T min, T max, T& value)
  {
    GetValue(name, value);
    if (value < min || value > max)
      throw string("Error in ExtStream::GetValue: the value of \"")
        + name + string("\" in \"") + file_name_ + "\" is "
        + to_str(value) + " but it should be in [" + to_str(min)
        + ", " + to_str(max) + "].";
  }

  /*! \brief Gets the value of a given variable without extracting them from
    the stream. */
  /*!
    Gets the (numerical) value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    Nothing is extracted from the stream.
    \param name the name of the variable.
    \param min the minimum value that the variable should take.
    \param max the maximum value that the variable should take.
    \param value value associated with the variable.
  */
  template <class T>
  void ExtStream::PeekValue(string name, T min, T max, T& value)
  {
    std::streampos initial_position = this->tellg();
    iostate state = this->rdstate();

    this->GetValue(name, min, max, value);

    this->clear(state);
    this->seekg(initial_position);
  }

  //! Gets the value of a given variable.
  /*!
    Gets the (numerical) value of a given variable, i.e. the next valid (not
    in a discarded line) number or element following the variable name. This
    methods also checks that the value meets given constraints.
    \param name the name of the variable.
    \param constraint the list of constraints. The constraints are delimited
    by |. The supported constraints are: positive, strictly positive,
    negative, strictly negative, non zero, integer, > x, >= x, < x, <= x, != x
    y z, = x y z.
    \param value value associated with the variable.
  */
  template <class T>
  void ExtStream::GetValue(string name, string constraint, T& value)
  {
    GetValue(name, value);
    if (!satisfies_constraint(value, constraint))
      throw string("Error in ExtStream::GetValue: the value of \"")
        + name + string("\" in \"") + file_name_ + "\" is "
        + to_str(value) + " but it should satisfy the following "
        + "constraint(s):\n" + show_constraint(constraint);
  }

  /*! \brief Gets the value of a given variable without extracting them from
    the stream. */
  /*!
    Gets the (numerical) value of a given variable, i.e. the next valid (not
    in a discarded line) number or element following the variable name. This
    methods also checks that the value meets given constraints. Nothing is
    extracted from the stream.
    \param name the name of the variable.
    \param constraint the list of constraints. The constraints are delimited
    by |. The supported constraints are: positive, strictly positive,
    negative, strictly negative, non zero, integer, > x, >= x, < x, <= x, != x
    y z, = x y z.
    \param value value associated with the variable.
  */
  template <class T>
  void ExtStream::PeekValue(string name, string constraint, T& value)
  {
    std::streampos initial_position = this->tellg();
    iostate state = this->rdstate();

    this->GetValue(name, constraint, value);

    this->clear(state);
    this->seekg(initial_position);
  }

  //! Gets the value of a given variable.
  /*!
    Gets the value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    \param name the name of the variable.
    \param value value associated with the variable.
  */
  void ExtStream::GetValue(string name, string& value)
  {
    SearchScope s(*this, name);

    string element;
    while (GetElement(element) && element != name);

    if (element != name)
      throw string("Error in ExtStream::GetValue: \"")
        + name + string("\" not found in \"") + file_name_ + "\".";

    if (!this->GetElement(value))
      throw string("Error in ExtStream::GetValue: ")
        + string("unable to get a value for \"") + name + string("\" in \"")
        + file_name_ + "\".";
  }

  /*! \brief Gets the value of a given variable without extracting them from
    the stream. */
  /*!
    Gets the value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    Nothing is extracted from the stream.
    \param name the name of the variable.
    \param value value associated with the variable.
  */
  void ExtStream::PeekValue(string name, string& value)
  {
    SearchScope s(*this, name);

    std::streampos initial_position = this->tellg();
    iostate state = this->rdstate();

    this->GetValue(name, value);

    this->clear(state);
    this->seekg(initial_position);
  }

  //! Gets the value of a given variable.
  /*!
    Gets the value of a given variable, i.e. the next valid (not in a
    discarded line) number or element following the variable name.  In
    addition, this method checks that the value is in an acceptable list of
    values.
    \param name the name of the variable.
    \param accepted list of accepted values.
    \param value value associated with the variable.
    \param delimiter delimiter in \a accepted. Default: |.
  */
  void ExtStream::GetValue(string name, string accepted, string& value,
                           string delimiter = "|")
  {
    GetValue(name, value);
    CheckAccepted(name, value, accepted, delimiter);
  }

  /*! \brief Gets the value of a given variable without extracting them from
    the stream. */
  /*!
    Gets the value of a given variable, i.e. the next valid (not in a
    discarded line) number or element following the variable name.  In
    addition, this method checks that the value is in an acceptable list of
    values. Nothing is extracted from the stream.
    \param name the name of the variable.
    \param accepted list of accepted values.
    \param value value associated with the variable.
    \param delimiter delimiter in \a accepted. Default: |.
  */
  void ExtStream::PeekValue(string name, string accepted, string& value,
                            string delimiter = "|")
  {
    SearchScope s(*this, name);

    std::streampos initial_position = this->tellg();
    iostate state = this->rdstate();

    this->GetValue(name, accepted, value, delimiter);

    this->clear(state);
    this->seekg(initial_position);
  }

  //! Gets the value of a given variable.
  /*!
    Gets the value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    \param name the name of the variable.
    \param value boolean associated with the variable.
  */
  void ExtStream::GetValue(string name, bool& value)
  {
    SearchScope s(*this, name);

    string element;
    while (GetElement(element) && element != name);

    if (element != name)
      throw string("Error in ExtStream::GetValue: \"")
        + name + string("\" not found in \"") + file_name_ + "\".";

    if (!this->GetElement(value))
      throw string("Error in ExtStream::GetValue: ")
        + string("unable to get a value for \"") + name + string("\" in \"")
        + file_name_ + "\".";
  }

  //! Gets the value of a given variable without extracting from the stream.
  /*!
    Gets the value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    Nothing is extracted from the stream.
    \param name the name of the variable.
    \param value boolean associated with the variable.
  */
  void ExtStream::PeekValue(string name, bool& value)
  {
    SearchScope s(*this, name);

    std::streampos initial_position = this->tellg();
    iostate state = this->rdstate();

    this->GetValue(name, value);

    this->clear(state);
    this->seekg(initial_position);
  }

  //! Checks that a value is in a given list of accepted values.
  /*!
    \param name the name of the entry with value \a value.
    \param value the value to be checked.
    \param accepted the list of accepted values.
    \param delimiter delimiter in \a accepted.
  */
  void ExtStream::CheckAccepted(string name, string value,
                                string accepted, string delimiter) const
  {
    vector<string> accepted_list = split(accepted, delimiter);
    int i = 0;
    while (i < int(accepted_list.size()) && trim(accepted_list[i]) != value)
      i++;
    if (i == int(accepted_list.size()))
      {
        string list = "[";
        for (i = 0; i < int(accepted_list.size()) - 1; i++)
          list += trim(accepted_list[i]) + " " + delimiter[0] + " ";
        if (accepted_list.size() != 0)
          list += trim(accepted_list[accepted_list.size() - 1]) + "]";
        throw string("Error in ExtStream::GetValue: the value of \"")
          + name + string("\" in \"") + file_name_ + "\" is \""
          + to_str(value) + "\" but it should be in " + list + ".";
      }
  }


  //////////////////
  // CONFIGSTREAM //
  //////////////////

  //! Default constructor.
  /*! Nothing is performed.
   */
  ConfigStream::ConfigStream(): ExtStream()
  {
    markup_tags_ = "<>$";
    section_ = "";
  }

  //! Main constructor.
  /*! Opens a file.
    \param file_name file to be opened.
  */
  ConfigStream::ConfigStream(string file_name,
                             string comments,
                             string delimiters,
                             string markup_tags):
    ExtStream(file_name, comments, delimiters),
    markup_tags_(markup_tags)
  {
    section_ = "";
  }

  //! Deselects the section.
  /*!
    Deselects the section (this is equivalent to SetSection("")).
  */
  void ConfigStream::NoSection()
  {
    section_ = "";
  }

  //! Sets the current section.
  /*!
    \param section current section.
  */
  void ConfigStream::SetSection(string section)
  {
    section_ = "";
    if (section != "")
      {
        this->FindFromBeginning(section);
        section_ = section;
      }
  }

  //! Returns the current section.
  /*!
    \return The current section.
  */
  string ConfigStream::GetSection() const
  {
    return section_;
  }

  //! Sets the markup tags.
  /*!
    \param markup_tags the new markup tags.
  */
  void ConfigStream::SetMarkupTags(string markup_tags)
  {
    markup_tags_ = markup_tags;
  }

  //! Returns the markup tags.
  /*!
    \return The markup tags.
  */
  string ConfigStream::GetMarkupTags() const
  {
    return markup_tags_;
  }

  //! Checks whether the stream or the current section is empty.
  /*!
    If no section has been selected, this method checks whether the stream has
    still valid elements to be read. If the stream is bound to a given
    section, this method checks whether there remains at least one element
    in the section.
    \return 'true' is the stream is empty, 'false' otherwise.
  */
  bool ConfigStream::IsEmpty()
  {
    std::streampos initial_position;
    bool success;

    initial_position = this->tellg();
    iostate state = this->rdstate();

    string element = ExtStream::GetElement();
    success = element != "";

    this->clear(state);
    this->seekg(initial_position);

    return !success || (section_ != "" && IsSection(element));
  }

  //! Checks whether the element is given.
  /*!
    \param element the element to be found.
    \return true if the element was found, false otherwise.
    \note The scope of the search is only the current section if any.
  */
  bool ConfigStream::Check(string element)
  {
    SearchScope s(*this, element);

    std::streampos initial_position = this->tellg();

    string elt;
    while (ExtStream::GetRawElement(elt) && elt != element
           && (section_ == "" || !IsSection(elt)));

    this->seekg(initial_position);
    return elt == element;
  }


  //! Sets the position of the get pointer after a given element.
  /*!
    Sets the position of the get pointer exactly after a given element.
    \param element the element to be found.
    \return true if the element was found, false otherwise.
    \note The scope of the search is only the current section if any.
  */
  bool ConfigStream::Find(string element)
  {
    SearchScope s(*this, element);

    string elt;
    while (ExtStream::GetRawElement(elt) && elt != element
           && (section_ == "" || !IsSection(elt)));

    if (section_ != "" && (elt == "" || IsSection(elt)))
      throw string("Error in ConfigStream::Find: end of section \"")
        + section_ + string("\" has been reached in file \"")
        + this->file_name_ + string("\".\nUnable to find \"")
        + element + "\".";
    if (elt == "")
      throw string("Error in ConfigStream::Find: \"")
        + element + string("\" not found in \"") + this->file_name_ + "\".";

    return elt == element;
  }

  //! Sets the position of the get pointer after a given element.
  /*!
    Sets the position of the get pointer exactly after a given element.
    \param element the element to be found from the beginning of the stream.
    \return true if the element was found, false otherwise.
    \note The current section (if any) is unset.
  */
  bool ConfigStream::FindFromBeginning(string element)
  {
    NoSection();
    this->Rewind();
    return this->Find(element);
  }

  //! Returns the next valid element.
  /*!
    Returns the next valid element, i.e. the next element that is
    not in a line to be discarded.
    \return The next valid element.
    \note Markups are replaced with their values.
  */
  string ConfigStream::GetElement()
  {
    string tmp;

    string element = ExtStream::GetElement();

    std::streampos initial_position = this->tellg();
    iostate state = this->rdstate();

    vector<string> elements;
    vector<bool> is_markup;

    split_markup(element, elements, is_markup, markup_tags_);

    element = "";

    for (int i = 0; i < int(elements.size()); i++)
      if (!is_markup[i])
        element += elements[i];
      else
        {
          this->Rewind();
          tmp = ExtStream::GetElement();
          while (tmp != elements[i] && tmp != "")
            tmp = ExtStream::GetElement();
          if (tmp == "")
            throw string("Error in ConfigStream::GetElement:")
              + string(" the value of the markup \"")
              + elements[i] + string("\" was not found in \"")
              + file_name_ + "\".";
          element += this->GetElement();
        }

    this->clear(state);
    this->seekg(initial_position);

    if (!section_.empty() && (element == "" || IsSection(element)))
      {
        string message = string("End of section \"") + section_
          + string("\" has been reached in file \"") + this->file_name_
          + "\".";
        if (this->searching_ != "")
          message += string("\nUnable to find \"")
            + this->searching_ + string("\".");
        throw message;
      }

    return element;
  }

  //! Returns the next valid line.
  /*!
    Returns the next valid line, i.e. the next line that is
    not a line to be discarded and from which comments have been extracted.
    \return The next valid line.
  */
  string ConfigStream::GetLine()
  {
    string tmp;

    string element = ExtStream::GetLine();

    std::streampos initial_position = this->tellg();
    iostate state = this->rdstate();

    vector<string> elements;
    vector<bool> is_markup;

    split_markup(element, elements, is_markup, markup_tags_);

    element = "";

    for (int i = 0; i < int(elements.size()); i++)
      if (!is_markup[i])
        element += elements[i];
      else
        {
          this->Rewind();
          tmp = ExtStream::GetElement();
          while (tmp != elements[i] && tmp != "")
            tmp = ExtStream::GetElement();
          if (tmp == "")
            throw string("Error in ConfigStream::GetLine:")
              + string(" the value of the markup \"")
              + elements[i] + string("\" was not found in \"")
              + file_name_ + "\".";
          element += this->GetElement();
        }

    this->clear(state);
    this->seekg(initial_position);

    return element;
  }

  //! Returns the next valid line.
  /*!
    Returns the next valid line, i.e. the next line that is
    not a line to be discarded and from which comments have been extracted.
    \param line (output) the next valid line.
  */
  bool ConfigStream::GetLine(string& line)
  {
    string tmp;

    bool success = ExtStream::GetLine(line);

    std::streampos initial_position = this->tellg();
    iostate state = this->rdstate();

    vector<string> elements;
    vector<bool> is_markup;

    split_markup(line, elements, is_markup, markup_tags_);

    line = "";

    for (int i = 0; i < int(elements.size()); i++)
      if (!is_markup[i])
        line += elements[i];
      else
        {
          this->Rewind();
          tmp = ExtStream::GetElement();
          while (tmp != elements[i] && tmp != "")
            tmp = ExtStream::GetElement();
          if (tmp == "")
            throw string("Error in ConfigStream::GetLine:")
              + string(" the value of the markup \"")
              + elements[i] + string("\" was not found in \"")
              + file_name_ + "\".";
          line += this->GetElement();
        }

    this->clear(state);
    this->seekg(initial_position);

    return success;
  }

  //! Checks whether a string is a section flag.
  /*!
    \param str string to be tested.
    \return True if 'str' is a section flag, false otherwise.
  */
  bool ConfigStream::IsSection(string str) const
  {
    return str[0] == '[' && str[str.size() - 1] == ']';
  }


  ///////////////////
  // CONFIGSTREAMS //
  ///////////////////

  //! Default constructor.
  /*! Nothing is performed.
   */
  ConfigStreams::ConfigStreams():
    streams_(0), current_(streams_.begin()), section_(""), searching_("")
  {
  }

  //! Main constructor.
  /*! Opens a set of file.
    \param files files to be opened.
  */
  ConfigStreams::ConfigStreams(const vector<string>& files):
    streams_(files.size()), current_(streams_.begin()), section_(""),
    searching_("")
  {
    for (int i = 0; i < int(files.size()); i++)
      streams_[i] = new ConfigStream(files[i]);
  }

  //! Constructor.
  /*! Opens a file.
    \param file file to be opened.
  */
  ConfigStreams::ConfigStreams(string file):
    streams_(1), current_(streams_.begin()), section_(""), searching_("")
  {
    streams_[0] = new ConfigStream(file);
  }

  //! Constructor.
  /*! Opens two files.
    \param file0 first file to be opened.
    \param file1 second file to be opened.
  */
  ConfigStreams::ConfigStreams(string file0, string file1):
    streams_(2), current_(streams_.begin()), section_(""), searching_("")
  {
    streams_[0] = new ConfigStream(file0);
    streams_[1] = new ConfigStream(file1);
  }

  //! Constructor.
  /*! Opens three files.
    \param file0 first file to be opened.
    \param file1 second file to be opened.
    \param file2 third file to be opened.
  */
  ConfigStreams::ConfigStreams(string file0, string file1, string file2):
    streams_(3), current_(streams_.begin()), section_(""), searching_("")
  {
    streams_[0] = new ConfigStream(file0);
    streams_[1] = new ConfigStream(file1);
    streams_[2] = new ConfigStream(file2);
  }

  //! Destructor.
  ConfigStreams::~ConfigStreams()
  {
    for (current_ = streams_.begin(); current_ != streams_.end(); ++current_)
      delete(*current_);
  }

  //! Returns the vector of streams.
  /*!
    \return A reference to the vector of ConfigStream.
  */
  vector<ConfigStream*>& ConfigStreams::GetStreams()
  {
    return streams_;
  }

  //! Returns the current position in the vector of streams.
  /*!
    \return An iterator that points to the current element of the
    vector of streams.
  */
  vector<ConfigStream*>::iterator ConfigStreams::GetCurrent()
  {
    return current_;
  }

  //! Adds a file to the streams.
  /*!
    \param file file to be added.
  */
  void ConfigStreams::AddFile(string file)
  {
    unsigned int l = current_ - streams_.begin();
    streams_.push_back(new ConfigStream(file));
    current_ = streams_.begin() + l;
  }

  //! Deselects the section.
  /*!
    Deselects the section (this is equivalent to SetSection("")) and rewinds
    the stream.
  */
  void ConfigStreams::NoSection()
  {
    section_ = "";
    for (current_ = streams_.begin(); current_ != streams_.end(); ++current_)
      (*current_)->SetSection(section_);
  }

  //! Sets the current section.
  /*!
    \param section current section.
  */
  void ConfigStreams::SetSection(string section)
  {
    section_ = section;
    for (current_ = streams_.begin(); current_ != streams_.end(); ++current_)
      (*current_)->SetSection("");

    current_ = streams_.begin();
    bool found = false;
    try
      {
        (*current_)->SetSection(section);
        found = true;
      }
    catch (...)
      {
        (*current_)->Rewind();
        (*current_)->section_ = section_;
      }
    while (!found && current_ != streams_.end() - 1)
      {
        ++current_;
        try
          {
            (*current_)->SetSection(section_);
            found = true;
          }
        catch (...)
          {
            (*current_)->Rewind();
            (*current_)->section_ = section_;
          }
      }

    if (!found)
      throw string("Error in ConfigStreams::SetSection: section \"")
        + section + string("\" was not found in ") + FileNames() + ".";
  }

  //! Returns the current section.
  /*!
    \return The current section.
  */
  string ConfigStreams::GetSection() const
  {
    return section_;
  }

  //! Checks whether a line should be discarded.
  /*!
    \param line line to be checked.
  */
  bool ConfigStreams::Discard(string line) const
  {
    return (*current_)->Discard(line);
  }

  //! Skips discarded lines.
  /*!
    Extracts discarded lines.
    \return A reference to 'this'.
  */
  ConfigStreams& ConfigStreams::SkipDiscarded()
  {
    (*current_)->SkipDiscarded();
    while (current_ != streams_.end() - 1 && !(*current_)->good())
      {
        ++current_;
        (*current_)->SkipDiscarded();
      }
    return *this;
  }

  //! Skips delimiters.
  /*!
    Extracts following delimiters from the string, until another character
    is found.
    \return A reference to this.
  */
  ConfigStreams& ConfigStreams::SkipDelimiters()
  {
    (*current_)->SkipDelimiters();
    while (current_ != streams_.end() - 1 && !(*current_)->good())
      {
        ++current_;
        (*current_)->SkipDelimiters();
      }
    return *this;
  }

  //! Removes delimiters at both ends of a string.
  /*!
    Removes delimiters at the beginning and at the end of a string.
    \param str string to processed.
    \return The string with delimiters removed.
  */
  string ConfigStreams::RemoveDelimiters(const string& str) const
  {
    return (*current_)->RemoveDelimiters(str);
  }

  //! Skips discarded lines and delimiters.
  /*!
    Extracts discarded lines and delimiters.
    \return A reference to this.
  */
  ConfigStreams& ConfigStreams::Skip()
  {
    this->SkipDiscarded();
    return this->SkipDelimiters();
  }

  //! Checks whether the streams are empty.
  /*!
    Checcks whether the streams have still valid elements to be read.
    \return 'true' is the streams are empty, 'false' otherwise.
  */
  bool ConfigStreams::IsEmpty()
  {
    return (*current_)->IsEmpty();
  }

  //! Rewinds all streams and goes back to the first stream.
  /*!
    Goes back to the beginning of each stream, clears the control state and
    goes back to the first stream.
    \return A reference to 'this'.
  */
  ConfigStreams& ConfigStreams::Rewind()
  {
    for (current_ = streams_.begin(); current_ != streams_.end(); ++current_)
      (*current_)->Rewind();
    current_ = streams_.begin();

    return *this;
  }

  //! Returns the next line.
  /*!
    \return The next line.
  */
  string ConfigStreams::GetFullLine()
  {
    string line;
    std::getline(**current_, line);

    if (is_emptystream(**current_) && current_ != streams_.end() - 1)
      ++current_;

    return line;
  }

  //! Returns the next line.
  /*!
    \param line (output) the next line.
  */
  bool ConfigStreams::GetFullLine(string& line)
  {
    bool success = bool(std::getline(**current_, line));

    if (is_emptystream(**current_) && current_ != streams_.end() - 1)
      ++current_;

    return success;
  }

  //! Returns the next line without extracting it from the stream.
  /*!
    \return The next line.
  */
  string ConfigStreams::PeekFullLine()
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    string line = (*current_)->PeekFullLine();

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);

    return line;
  }

  //! Returns the next line without extracting it from the stream.
  /*!
    \param position (output) the position of the line following the next line.
    \return The next line.
  */
  string ConfigStreams::PeekFullLine(std::streampos& position)
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    string line = (*current_)->PeekFullLine(position);

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);

    return line;
  }

  //! Returns the next line without extracting it from the stream.
  /*!
    \param line (output) the next line.
  */
  bool ConfigStreams::PeekFullLine(string& line)
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    bool success = (*current_)->PeekFullLine(line);

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);

    return success;
  }

  //! Skips full lines.
  /*!
    \param nb number of lines to be skipped.
  */
  void ConfigStreams::SkipFullLines(int nb)
  {
    for (int i = 0; i < nb; i++)
      this->GetFullLine();
  }

  //! Returns the next valid line, without markups substitution.
  /*!
    Returns the next valid line, i.e. the next line that is
    not a line to be discarded and from which comments have been extracted.
    No markup substitution is performed.
    \return The next valid line.
  */
  string ConfigStreams::GetRawLine()
  {
    this->Skip();
    return (*current_)->ExtStream::GetLine();
  }

  //! Returns the next valid line.
  /*!
    Returns the next valid line, i.e. the next line that is
    not a line to be discarded and from which comments have been extracted.
    \return The next valid line.
  */
  string ConfigStreams::GetLine()
  {
    string line;
    this->GetLine(line);
    return line;
  }

  //! Returns the next valid line.
  /*!
    Returns the next valid line, i.e. the next line that is
    not a line to be discarded and from which comments have been extracted.
    \param line (output) the next valid line.
  */
  bool ConfigStreams::GetLine(string& line)
  {
    string tmp;

    line = this->GetRawLine();
    bool success = (line != "");

    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    vector<string> elements;
    vector<bool> is_markup;

    split_markup(line, elements, is_markup, (*current_)->GetMarkupTags());

    line = "";

    for (int i = 0; i < int(elements.size()); i++)
      if (!is_markup[i])
        line += elements[i];
      else
        {
          this->Rewind();
          tmp = this->GetRawElement();
          while (tmp != elements[i] && tmp != "")
            tmp = this->GetRawElement();
          if (tmp == "")
            throw string("Error in ConfigStreams::GetLine:")
              + string(" the value of the markup \"")
              + elements[i] + string("\" was not found in ")
              + FileNames() + ".";
          line += this->GetElement();
        }

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);

    return success;
  }

  //! Returns the next valid line without extracting it from the stream.
  /*!
    Returns the next valid line, i.e. the next line that is
    not a line to be discarded and from which comments have been extracted.
    Nothing is extracted from the stream.
    \return The next valid line.
  */
  string ConfigStreams::PeekLine()
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    string line = this->GetLine();

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);

    return line;
  }

  //! Returns the next valid line without extracting it from the stream.
  /*!
    Returns the next valid line, i.e. the next line that is
    not a line to be discarded and from which comments have been extracted.
    Nothing is extracted from the stream.
    \param position (output) the position of the line following the next valid
    line.
    \return The valid line.
  */
  string ConfigStreams::PeekLine(std::streampos& position)
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    string line = this->GetLine();
    position = (*current_)->tellg();

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);

    return line;
  }

  //! Returns the next valid line without extracting it from the stream.
  /*!
    Returns the next valid line, i.e. the next line that is
    not a line to be discarded and from which comments have been extracted.
    Nothing is extracted from the stream.
    \param line (output) the next valid line.
  */
  bool ConfigStreams::PeekLine(string& line)
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    bool success = this->GetLine(line);

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);

    return success;
  }

  //! Skips valid lines.
  /*!
    \param nb number of lines to be skipped.
  */
  void ConfigStreams::SkipLines(int nb)
  {
    for (int i = 0; i < nb; i++)
      this->GetLine();
  }

  //! Sets the position of the get pointer after a given element.
  /*!
    Sets the position of the get pointer exactly after a given element.
    \param element the element to be found.
    \return true if the element was found, false otherwise.
  */
  bool ConfigStreams::Find(string element)
  {
    SearchScope s(*this, element);

    bool found;
    try
      {
        found = (*current_)->Find(element);
      }
    catch (...)
      {
        found = false;
      }
    while (!found && current_ != streams_.end() - 1)
      {
        ++current_;
        try
          {
            found = (*current_)->Find(element);
          }
        catch (...)
          {
            found = false;
          }
      }
    if (!found && !section_.empty())
      throw string("Error in ConfigStreams::Find: end of section \"")
        + section_ + string("\" has been reached in ")
        + FileNames() + string(".\nUnable to find \"")
        + element + "\".";
    if (!found)
      throw string("Error in ConfigStreams::Find: \"")
        + element + string("\" not found in ") + FileNames() + ".";

    return found;
  }

  //! Sets the position of the get pointer after a given element.
  /*!
    Sets the position of the get pointer exactly after a given element.
    \param element the element to be found from the beginning of the stream.
    \return true if the element was found, false otherwise.
  */
  bool ConfigStreams::FindFromBeginning(string element)
  {
    NoSection();
    this->Rewind();
    return this->Find(element);
  }

  //! Returns the next valid element, without markups substitution.
  /*!
    Returns the next valid element, i.e. the next element that is
    not in a line to be discarded. No markup substitution is performed.
    \return The next valid element.
  */
  string ConfigStreams::GetRawElement()
  {
    string element;
    while ((element = (*current_)->ExtStream::GetElement()) == ""
           && current_ != streams_.end() - 1)
      ++current_;
    return element;
  }

  //! Returns the next valid element.
  /*!
    Returns the next valid element, i.e. the next element that is
    not in a line to be discarded.
    \return The next valid element.
  */
  string ConfigStreams::GetElement()
  {
    string tmp;

    string element = GetRawElement();

    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    vector<string> elements;
    vector<bool> is_markup;

    split_markup(element, elements, is_markup, (*current_)->GetMarkupTags());

    element = "";

    for (int i = 0; i < int(elements.size()); i++)
      if (!is_markup[i])
        element += elements[i];
      else
        {
          this->Rewind();
          tmp = GetRawElement();
          while (tmp != elements[i] && tmp != "")
            tmp = GetRawElement();
          if (tmp == "")
            throw string("Error in ConfigStreams::GetElement: ")
              + string("the value of the markup \"")
              + elements[i] + string("\" was not found in ")
              + FileNames() + ".";
          element += this->GetElement();
        }

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);

    if (!section_.empty() && (element == "" || IsSection(element)))
      {
        string message = string("End of section \"") + section_
          + string("\" has been reached in ") + FileNames() + ".";
        if (searching_ != "")
          message += string("\nUnable to find \"")
            + this->searching_ + string("\".");
        throw message;
      }

    return element;
  }

  //! Gets the next valid element.
  /*!
    Gets the next valid element, i.e. the next element that is
    not in a line to be discarded.
    \param element (output) the next valid element.
  */
  template <class T>
  bool ConfigStreams::GetElement(T& element)
  {
    string str = this->GetElement();
    convert(str, element);

    return (str != "");
  }

  //! Returns the next valid element without extracting it from the stream.
  /*!
    Returns the next valid element, i.e. the next element that is
    not in a line to be discarded.
    \return The next valid element.
  */
  string ConfigStreams::PeekElement()
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    string element = this->GetElement();

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);

    return element;
  }

  //! Gets the next valid element without extracting it from the stream.
  /*!
    Gets the next valid element, i.e. the next element that is
    not in a line to be discarded.
    \param element (output) the next valid element.
  */
  template <class T>
  bool ConfigStreams::PeekElement(T& element)
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    bool success = this->GetElement(element);

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);

    return success;
  }

  //! Skips valid elements.
  /*!
    \param nb number of valid elements to be skipped.
  */
  void ConfigStreams::SkipElements(int nb)
  {
    for (int i = 0; i < nb; i++)
      this->GetElement();
  }

  //! Returns the next valid number.
  /*!
    Returns the next valid number, i.e. the next number that is
    not in a line to be discarded.
    \return The next valid number.
  */
  double ConfigStreams::GetNumber()
  {
    string element;
    while (this->GetElement(element) && !is_num(element));

    return is_num(element) ? to_num<double>(element) : 0.;
  }

  //! Gets the next valid number.
  /*!
    Gets the next valid number, i.e. the next number that is
    not in a line to be discarded.
    \param element (output) the next valid number.
  */
  template <class T>
  bool ConfigStreams::GetNumber(T& number)
  {
    string element;
    bool success;
    while ((success = this->GetElement(element)) && !is_num(element));

    number = is_num(element) ? to_num<T>(element) : T(0);

    return success;
  }

  //! Returns the next valid number without extracting it from the stream.
  /*!
    Returns the next valid number, i.e. the next number that is
    not in a line to be discarded.
    \return The next valid number.
  */
  double ConfigStreams::PeekNumber()
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    double num = this->GetNumber();

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);

    return num;
  }

  //! Gets the next valid number without extracting it from the stream.
  /*!
    Gets the next valid number, i.e. the next number that is
    not in a line to be discarded.
    \param number (output) the next valid number.
  */
  template <class T>
  bool ConfigStreams::PeekNumber(T& number)
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    bool success = this->GetNumber(number);

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);

    return success;
  }

  //! Skips numbers.
  /*!
    \param nb number of numbers to be skipped.
  */
  void ConfigStreams::SkipNumbers(int nb)
  {
    for (int i = 0; i < nb; i++)
      this->GetNumber();
  }

  //! Gets the value of a given variable.
  /*!
    Gets the value of a given variable, i.e. the next valid
    (not in a discarded line) element following the variable name.
    \param name the name of the variable.
    \return the value of the variable.
  */
  string ConfigStreams::GetValue(string name)
  {
    SearchScope s(*this, name);

    string element;
    while (this->GetElement(element) && element != name);

    if (element != name)
      throw string("Error in ConfigStreams::GetValue: \"")
        + name + string("\" not found in ") + FileNames() + ".";

    return this->GetElement();
  }

  //! Gets the value of a given variable without extracting from the stream.
  /*!
    Gets the value of a given variable, i.e. the next valid
    (not in a discarded line) element following the variable name.
    Nothing is extracted from the stream.
    \param name the name of the variable.
    \return the value associated with the variable.
  */
  string ConfigStreams::PeekValue(string name)
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    string element = this->GetValue(name);

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);

    return element;
  }

  //! Gets the value of a given variable.
  /*!
    Gets the (numerical) value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    \param name the name of the variable.
    \param value value associated with the variable.
  */
  template <class T>
  void ConfigStreams::GetValue(string name, T& value)
  {
    SearchScope s(*this, name);

    string element;
    while (this->GetElement(element) && element != name);

    if (element != name)
      throw string("Error in ConfigStreams::GetValue: \"")
        + name + string("\" not found in ") + FileNames() + ".";
    if (!this->GetElement(element))
      throw string("Error in ConfigStreams::GetValue: unable to read value")
        + string(" of \"") + name + string("\" in ") + FileNames() + ".";
    if (!is_num(element))
      throw string("Error in ConfigStreams::GetValue: the value of \"") + name
        + string("\" in ") + FileNames() + string(" is \"") + element
        + "\", but it should be a number.";

    value = to_num<T>(element);
  }

  //! Gets the value of a given variable.
  /*!
    Gets the integral value of a given variable, i.e. the next valid (not in a
    discarded line) number or element following the variable name.
    \param name the name of the variable.
    \param value value associated with the variable.
  */
  void ConfigStreams::GetValue(string name, int& value)
  {
    SearchScope s(*this, name);

    string element;
    while (this->GetElement(element) && element != name);

    if (element != name)
      throw string("Error in ConfigStreams::GetValue: \"")
        + name + string("\" not found in ") + FileNames() + ".";
    if (!this->GetElement(element))
      throw string("Error in ConfigStreams::GetValue: unable to read value")
        + string(" of \"") + name + string("\" in ") + FileNames() + ".";
    if (!is_integer(element))
      throw string("Error in ConfigStreams::GetValue: the value of \"") + name
        + string("\" in ") + FileNames() + string(" is \"") + element
        + "\", but it should be an integer.";

    value = to_num<int>(element);
  }

  /*! \brief Gets the value of a given variable without extracting them from
    the stream. */
  /*!
    Gets the (numerical) value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    Nothing is extracted from the stream.
    \param name the name of the variable.
    \param value value associated with the variable.
  */
  template <class T>
  void ConfigStreams::PeekValue(string name, T& value)
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    this->GetValue(name, value);

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);
  }

  //! Gets the value of a given variable.
  /*!
    Gets the (numerical) value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    \param name the name of the variable.
    \param min the minimum value that the variable should take.
    \param max the maximum value that the variable should take.
    \param value value associated with the variable.
  */
  template <class T>
  void ConfigStreams::GetValue(string name, T min, T max, T& value)
  {
    GetValue(name, value);
    if (value < min || value > max)
      throw string("Error in ConfigStreams::GetValue: the value of \"")
        + name + string("\" in ") +  FileNames() + " is "
        + to_str(value) + " but it should be in [" + to_str(min)
        + ", " + to_str(max) + "].";
  }

  /*! \brief Gets the value of a given variable without extracting them from
    the stream. */
  /*!
    Gets the (numerical) value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    Nothing is extracted from the stream.
    \param name the name of the variable.
    \param min the minimum value that the variable should take.
    \param max the maximum value that the variable should take.
    \param value value associated with the variable.
  */
  template <class T>
  void ConfigStreams::PeekValue(string name, T min, T max, T& value)
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    this->GetValue(name, min, max, value);

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);
  }

  //! Gets the value of a given variable.
  /*!
    Gets the (numerical) value of a given variable, i.e. the next valid (not
    in a discarded line) number or element following the variable name. This
    methods also checks that the value meets given constraints.
    \param name the name of the variable.
    \param constraint the list of constraints. The constraints are delimited
    by |. The supported constraints are: positive, strictly positive,
    negative, strictly negative, non zero, integer, > x, >= x, < x, <= x, != x
    y z, = x y z.
    \param value value associated with the variable.
  */
  template <class T>
  void ConfigStreams::GetValue(string name, string constraint, T& value)
  {
    GetValue(name, value);
    if (!satisfies_constraint(value, constraint))
      throw string("Error in ConfigStreams::GetValue: the value of \"")
        + name + string("\" in ") + FileNames() + " is "
        + to_str(value) + " but it should satisfy the following "
        + "constraint(s):\n" + show_constraint(constraint);
  }

  /*! \brief Gets the value of a given variable without extracting them from
    the stream. */
  /*!
    Gets the (numerical) value of a given variable, i.e. the next valid (not
    in a discarded line) number or element following the variable name. This
    methods also checks that the value meets given constraints. Nothing is
    extracted from the stream.
    \param name the name of the variable.
    \param constraint the list of constraints. The constraints are delimited
    by |. The supported constraints are: positive, strictly positive,
    negative, strictly negative, non zero, integer, > x, >= x, < x, <= x, != x
    y z, = x y z.
    \param value value associated with the variable.
  */
  template <class T>
  void ConfigStreams::PeekValue(string name, string constraint, T& value)
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    this->GetValue(name, constraint, value);

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);
  }

  //! Gets the value of a given variable.
  /*!
    Gets the value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    \param name the name of the variable.
    \param value value associated with the variable.
  */
  void ConfigStreams::GetValue(string name, string& value)
  {
    SearchScope s(*this, name);

    string element;
    while (this->GetElement(element) && element != name);

    if (element != name)
      throw string("Error in ConfigStreams::GetValue: \"")
        + name + string("\" not found in ") + FileNames() + ".";

    if (!this->GetElement(value))
      throw string("Error in ConfigStreams::GetValue: ")
        + string("unable to get a value for \"") + name + string("\" in ")
        + FileNames() + ".";
  }

  /*! \brief Gets the value of a given variable without extracting them from
    the stream. */
  /*!
    Gets the value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    Nothing is extracted from the stream.
    \param name the name of the variable.
    \param value value associated with the variable.
  */
  void ConfigStreams::PeekValue(string name, string& value)
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    this->GetValue(name, value);

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);
  }

  //! Gets the value of a given variable.
  /*!
    Gets the value of a given variable, i.e. the next valid (not in a
    discarded line) number or element following the variable name.  In
    addition, this method checks that the value is in an acceptable list of
    values.
    \param name the name of the variable.
    \param accepted list of accepted values.
    \param value value associated with the variable.
    \param delimiter delimiter in \a accepted. Default: |.
  */
  void ConfigStreams::GetValue(string name, string accepted, string& value,
                               string delimiter = "|")
  {
    GetValue(name, value);
    CheckAccepted(name, value, accepted, delimiter);
  }

  /*! \brief Gets the value of a given variable without extracting them from
    the stream. */
  /*!
    Gets the value of a given variable, i.e. the next valid (not in a
    discarded line) number or element following the variable name.  In
    addition, this method checks that the value is in an acceptable list of
    values. Nothing is extracted from the stream.
    \param name the name of the variable.
    \param accepted list of accepted values.
    \param value value associated with the variable.
    \param delimiter delimiter in \a accepted. Default: |.
  */
  void ConfigStreams::PeekValue(string name, string accepted, string& value,
                                string delimiter = "|")
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    this->GetValue(name, accepted, value, delimiter);

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);
  }

  //! Gets the value of a given variable.
  /*!
    Gets the value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    \param name the name of the variable.
    \param value boolean associated with the variable.
  */
  void ConfigStreams::GetValue(string name, bool& value)
  {
    SearchScope s(*this, name);

    string element;
    while (this->GetElement(element) && element != name);

    if (element != name)
      throw string("Error in ConfigStreams::GetValue: \"")
        + name + string("\" not found in ") + FileNames() + ".";

    if (!this->GetElement(value))
      throw string("Error in ConfigStreams::GetValue: ")
        + string("unable to get a value for \"") + name + string("\" in ")
        + FileNames() + ".";
  }

  //! Gets the value of a given variable without extracting from the stream.
  /*!
    Gets the value of a given variable, i.e. the next valid
    (not in a discarded line) number or element following the variable name.
    Nothing is extracted from the stream.
    \param name the name of the variable.
    \param value boolean associated with the variable.
  */
  void ConfigStreams::PeekValue(string name, bool& value)
  {
    vector<ConfigStream*>::iterator iter = current_;
    std::streampos initial_position = (*current_)->tellg();
    ifstream::iostate state = (*current_)->rdstate();

    this->GetValue(name, value);

    this->Rewind();
    current_ = iter;
    (*current_)->clear(state);
    (*current_)->seekg(initial_position);
  }

  //! Checks whether a string is a section flag.
  /*!
    \param str string to be tested.
    \return True if 'str' is a section flag, false otherwise.
  */
  bool ConfigStreams::IsSection(string str) const
  {
    return str[0] == '[' && str[str.size() - 1] == ']';
  }

  //! Returns file names in string form.
  /*!
    \return File names in order to print them on screen.
  */
  string ConfigStreams::FileNames() const
  {
    string output = "";
    int Nstream = streams_.size();
    for (int i = 0; i < Nstream - 2; i++)
      output += string("\"") + streams_[i]->GetFileName() + "\", ";
    if (Nstream > 1)
      output += string("\"") + streams_[Nstream - 2]->GetFileName()
        + "\" or ";
    if (Nstream > 0)
      output += string("\"") + streams_[Nstream - 1]->GetFileName() + "\"";
    return output;
  }

  //! Checks that a value is in a given list of accepted values.
  /*!
    \param name the name of the entry with value \a value.
    \param value the value to be checked.
    \param accepted the list of accepted values.
    \param delimiter delimiter in \a accepted.
  */
  void ConfigStreams::CheckAccepted(string name, string value,
                                    string accepted, string delimiter) const
  {
    vector<string> accepted_list = split(accepted, delimiter);
    int i = 0;
    while (i < int(accepted_list.size()) && trim(accepted_list[i]) != value)
      i++;
    if (i == int(accepted_list.size()))
      {
        string list = "[";
        for (i = 0; i < int(accepted_list.size()) - 1; i++)
          list += trim(accepted_list[i]) + " " + delimiter[0] + " ";
        if (accepted_list.size() != 0)
          list += trim(accepted_list[accepted_list.size() - 1]) + "]";
        throw string("Error in ConfigStreams::GetValue: the value of \"")
          + name + string("\" in ") + FileNames() + " is \""
          + to_str(value) + "\" but it should be in " + list + ".";
      }
  }

}  // namespace Talos.


#define TALOS_FILE_FILES_CXX
#endif
