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


#ifndef TALOS_FILE_DATE_CXX

#include "String.hxx"
#include "Date.hxx"

namespace Talos
{

  //////////
  // DATE //
  //////////

  //! Default constructor.
  Date::Date():
    year_(0), month_(0), day_(0),
    hour_(0), minutes_(0), seconds_(0),
    month_lengths_(12)
  {
    int month_lengths[12] = {31, 28, 31, 30, 31, 30,
                             31, 31, 30, 31, 30, 31
    };
    for (int i = 0; i < 12; ++i)
      month_lengths_[i] = month_lengths[i];

    this->Adjust();
  }

  //! Copy constructor.
  /*!
    \param date date to be copied.
  */
  Date::Date(const Date& date):
    year_(date.GetYear()), month_(date.GetMonth()), day_(date.GetDay()),
    hour_(date.GetHour()), minutes_(date.GetMinutes()),
    seconds_(date.GetSeconds()), month_lengths_(12)
  {
    int month_lengths[12] = {31, 28, 31, 30, 31, 30,
                             31, 31, 30, 31, 30, 31
    };
    for (int i = 0; i < 12; ++i)
      month_lengths_[i] = month_lengths[i];

    this->Adjust();
  }

  //! Copy constructor.
  /*!
    \param date  date in format YYYY, YYYYMM, YYYYMMDD, YYYYMMDDHH,
    YYYYMMDDHHMM or YYYYMMDDHHMMSS, with each item (YYYY, MM, DD, etc.)
    possibly delimited with any character(s).
  */
  Date::Date(string date): month_lengths_(12)
  {
    int month_lengths[12] = {31, 28, 31, 30, 31, 30,
                             31, 31, 30, 31, 30, 31
    };
    for (int i = 0; i < 12; ++i)
      month_lengths_[i] = month_lengths[i];

    this->SetDate(date);
  }

  //! Constructor.
  /*!
    \param yyyymmdd date in format YYYYMMDD.
  */
  Date::Date(int yyyymmdd):
    hour_(0), minutes_(0), seconds_(0),
    month_lengths_(12)
  {
    year_ = yyyymmdd / 10000;
    month_ = (yyyymmdd % 10000) / 100;
    day_ = yyyymmdd % 100;

    int month_lengths[12] = {31, 28, 31, 30, 31, 30,
                             31, 31, 30, 31, 30, 31
    };
    for (int i = 0; i < 12; ++i)
      month_lengths_[i] = month_lengths[i];

    this->Adjust();
  }

  //! Constructor.
  /*!
    \param yyyy year.
    \param mm month.
    \param dd day.
    \param hh hour.
    \param mn minutes.
    \param sc seconds.
  */
  Date::Date(int yyyy, int mm, int dd,
             int hh, int mn, double sc):
    year_(yyyy), month_(mm), day_(dd),
    hour_(hh), minutes_(mn), seconds_(sc),
    month_lengths_(12)
  {
    int month_lengths[12] = {31, 28, 31, 30, 31, 30,
                             31, 31, 30, 31, 30, 31
    };
    for (int i = 0; i < 12; ++i)
      month_lengths_[i] = month_lengths[i];

    this->Adjust();
  }

  //! Assignment operator.
  /*!
    \param date date to be copied.
  */
  Date& Date::operator=(const Date& date)
  {
    year_ = date.GetYear();
    month_ = date.GetMonth();
    day_ = date.GetDay();
    hour_ = date.GetHour();
    minutes_ = date.GetMinutes();
    seconds_ = date.GetSeconds();

    this->Adjust();

    return *this;
  }

  //! Assignment operator.
  /*!
    \param date  date in format YYYY, YYYYMM, YYYYMMDD, YYYYMMDDHH,
    YYYYMMDDHHMM or YYYYMMDDHHMMSS, with each item (YYYY, MM, DD, etc.)
    possibly delimited with any character(s).
  */
  Date& Date::operator=(string date)
  {
    SetDate(date);

    return *this;
  }

  //! Sets the date.
  /*!
    \param yyyymmdd date in format YYYYMMDD.
  */
  void Date::SetDate(int yyyymmdd)
  {
    year_ = yyyymmdd / 10000;
    month_ = (yyyymmdd % 10000) / 100;
    day_ = yyyymmdd % 100;

    hour_ = 0;
    minutes_ = 0;
    seconds_ = 0.;

    this->Adjust();
  }

  //! Sets the date.
  /*!
    \param date  date in format YYYY, YYYYMM, YYYYMMDD, YYYYMMDDHH,
    YYYYMMDDHHMM or YYYYMMDDHHMMSS, with each item (YYYY, MM, DD, etc.)
    possibly delimited with any character(s).
  */
  void Date::SetDate(string date)
  {
    vector<string> split_date;
    string item = "";
    int code;
    // Parses the string to split items.
    for (unsigned int i = 0; i < date.size(); i++)
      {
        code = int(date[i]);
        // If the current character is an integer.
        if (code > 47 && code < 58)
          item += date[i];
        else if (!item.empty())
          {
            split_date.push_back(item);
            item = "";
          }
      }
    if (!item.empty())
      split_date.push_back(item);

    // Checks that the date can be parsed.
    string compressed_date;
    if (split_date.size() == 0)
      throw string("Badly formatted date: \"") + date + string("\".");
    if (split_date[0].size() < 4)
      throw string("Badly formatted date: \"") + date + string("\".");
    for (unsigned int i = 0; i < split_date.size(); i++)
      {
        if (split_date[i].size() % 2 != 0 || !is_integer(split_date[i]))
          throw string("Badly formatted date: \"") + date + string("\".");
        compressed_date += split_date[i];
      }

    // Year is at least provided.
    year_ = to_num<int>(compressed_date.substr(0, 4));

    // Initialization (in case the date is not specified up to seconds).
    month_ = 1;
    day_ = 1;
    hour_ = 0;
    minutes_ = 0;
    seconds_ = 0.;

    // Retrieves all information.
    unsigned int length = compressed_date.size();
    if (length >= 6)
      month_ = to_num<int>(compressed_date.substr(4, 2));
    if (length >= 8)
      day_ = to_num<int>(compressed_date.substr(6, 2));
    if (length >= 10)
      hour_ = to_num<int>(compressed_date.substr(8, 2));
    if (length >= 12)
      minutes_ = to_num<int>(compressed_date.substr(10, 2));
    if (length >= 14)
      seconds_ = to_num<double>(compressed_date.substr(12, 2));

    if (!IsValid())
      throw string("Date \"") + date + string("\" is invalid.");

    this->Adjust();
  }

  //! Sets the date.
  /*!
    \param yyyy year.
    \param mm month.
    \param dd day.
    \param hh hour.
    \param mn minutes.
    \param sc seconds.
  */
  void Date::SetDate(int yyyy, int mm, int dd,
                     int hh, int mn, double sc)
  {
    year_ = yyyy;
    month_ = mm;
    day_ = dd;
    hour_ = hh;
    minutes_ = mn;
    seconds_ = sc;

    this->Adjust();
  }

  //! Adjusts month lengths according to the year.
  void Date::LeapYearAdjust()
  {
    if (LeapYear())
      month_lengths_[1] = 29;
    else
      month_lengths_[1] = 28;
  }

  //! Checks whether a date is valid.
  /*!
    \return True is the current date is valid, false otherwise.
  */
  bool Date::IsValid()
  {
    this->LeapYearAdjust();

    // Month.
    return month_ > 0 && month_ < 13
                  && day_ > 0 && day_ < month_lengths_[month_ - 1] + 1
                            && minutes_ > -1 && minutes_ < 60
                                               && seconds_ >= 0. && seconds_ < 60.;
  }

  //! Adjusts the date to make it valid.
  void Date::Adjust()
  {
    // Minutes.
    if (seconds_ >= 60.)
      {
        minutes_ += int(seconds_ / 60.);
        seconds_ -= double(int(seconds_ / 60.)) * 60.;
      }
    if (seconds_ < 0.)
      {
        minutes_ += int(seconds_ / 60.) - 1;
        seconds_ -= double(int(seconds_ / 60.) - 1) * 60.;
        if (seconds_ == 60.)
          {
            minutes_ += 1;
            seconds_ = 0.;
          }
      }

    // Hours.
    if (minutes_ > 59)
      {
        hour_ += minutes_ / 60;
        minutes_ = minutes_ % 60;
      }
    else if (minutes_ < 0)
      {
        hour_ += (minutes_ + 1) / 60 - 1;
        minutes_ = 59 + (minutes_ + 1) % 60;
      }

    // Days.
    if (hour_ > 23)
      {
        day_ += hour_ / 24;
        hour_ = hour_ % 24;
      }
    else if (hour_ < 0)
      {
        day_ += hour_ / 24 - 1;
        hour_ = 24 + hour_ % 24;
      }

    // Months.
    if (month_ > 12)
      {
        year_ += month_ / 12;
        month_ = month_ % 12;
      }
    else if (month_ <= 0)
      {
        year_ += month_ / 12 - 1;
        month_ = 12 + month_ % 12;
      }
    this->LeapYearAdjust();

    // Months again.
    while (day_ > month_lengths_[month_ - 1])
      {
        day_ -= month_lengths_[month_ - 1];
        ++month_;
        if (month_ == 13)
          {
            month_ = 1;
            ++year_;
            this->LeapYearAdjust();
          }
      }

    while (day_ <= 0)
      {
        --month_;
        if (month_ == 0)
          {
            month_ = 12;
            --year_;
            this->LeapYearAdjust();
          }
        day_ += month_lengths_[month_ - 1];
      }

  }

  //! Is a given year a leap year?
  /*!
    \param year year.
    \return true if the year 'year' is a leap year, false otherwise.
  */
  bool Date::LeapYear(int year) const
  {
    return (year % 4 == 0
            && (year % 100 != 0 || year % 400 == 0));
  }

  //! Is the current year a leap year?
  /*!
    \return true if the current year is a leap year, false otherwise.
  */
  bool Date::LeapYear() const
  {
    return this->LeapYear(year_);
  }

  //! Returns the date in format YYYYMMDD.
  /*!
    \return The date in format YYYYMMDD.
  */
  int Date::GetDate() const
  {
    return year_ * 10000 + month_ * 100 + day_;
  }

  //! Returns the date in a given format.
  /*!
    'format' defines the format of the output. Special sequences are
    %y, %m, %d, %h, %i, %s for the year, the month, the day, the hour,
    the minutes and the seconds respectively.
    \param format format.
    \return The date in format 'format'.
  */
  string Date::GetDate(const string& format) const
  {
    string output("");
    string::size_type index_b(0), index_e, tmp;

    while ((index_b != string::npos)
           && (index_b < format.size()))
      {
        tmp = format.substr(index_b).find_first_of("%");
        index_e = tmp == string::npos ? string::npos : index_b + tmp;
        if ((index_e != string::npos) && (index_e != format.size() - 1))
          {
            output += format.substr(index_b, index_e - index_b);
            if (format[index_e + 1] == 'y')
              output += to_str_fill(year_, 4, '0', ostringstream::right);
            else if (format[index_e + 1] == 'm')
              output += to_str_fill(month_, 2, '0', ostringstream::right);
            else if (format[index_e + 1] == 'd')
              output += to_str_fill(day_, 2, '0', ostringstream::right);
            else if (format[index_e + 1] == 'h')
              output += to_str_fill(hour_, 2, '0', ostringstream::right);
            else if (format[index_e + 1] == 'i')
              output += to_str_fill(minutes_, 2, '0', ostringstream::right);
            else if (format[index_e + 1] == 's')
              output += to_str_fill(seconds_, 2, '0', ostringstream::right);
            else
              output += format.substr(index_e, 2);
            index_b = index_e + 2;
          }
        else
          {
            output += format.substr(index_b);
            index_b = string::npos;
          }
      }

    return output;
  }

  //! Returns the date in a given format.
  /*!
    'format' defines the format of the output. Special sequences are
    %y, %m, %d, %h, %i, %s for the year, the month, the day, the hour,
    the minutes and the seconds respectively.
    \param format format.
    \param date (output) the date in format 'format'.
  */
  template <class T>
  void Date::GetDate(const string& format, T& date) const
  {
    date = to_num<T>(this->GetDate(format));
  }

  //! Returns the year.
  /*!
    \return The year.
  */
  int Date::GetYear() const
  {
    return year_;
  }

  //! Returns the month.
  /*!
    \return The month.
  */
  int Date::GetMonth() const
  {
    return month_;
  }

  //! Returns the day.
  /*!
    \return The day.
  */
  int Date::GetDay() const
  {
    return day_;
  }

  //! Returns the hour.
  /*!
    \return The hour.
  */
  int Date::GetHour() const
  {
    return hour_;
  }

  //! Returns the minutes.
  /*!
    \return The minutes.
  */
  int Date::GetMinutes() const
  {
    return minutes_;
  }

  //! Returns the seconds.
  /*!
    \return The seconds.
  */
  double Date::GetSeconds() const
  {
    return seconds_;
  }

  //! Adds years to the current date.
  /*!
    \param nb_yy number of years.
  */
  void Date::AddYears(int nb_yy)
  {
    year_ += nb_yy;
    this->Adjust();
  }

  //! Adds months to the current date.
  /*!
    \param nb_mm number of months.
  */
  void Date::AddMonths(int nb_mm)
  {
    month_ += nb_mm;
    this->Adjust();
  }

  //! Adds days to the current date.
  /*!
    \param nb_dd number of days.
  */
  void Date::AddDays(int nb_dd)
  {
    day_ += nb_dd;
    this->Adjust();
  }

  //! Adds hours to the current date.
  /*!
    \param nb_hh number of hours.
  */
  void Date::AddHours(int nb_hh)
  {
    hour_ += nb_hh;
    this->Adjust();
  }

  //! Adds minutes to the current date.
  /*!
    \param nb_mn number of minutes.
  */
  void Date::AddMinutes(int nb_mn)
  {
    minutes_ += nb_mn;
    this->Adjust();
  }

  //! Adds seconds to the current date.
  /*!
    \param nb_sc number of seconds.
  */
  void Date::AddSeconds(double nb_sc)
  {
    seconds_ += nb_sc;
    this->Adjust();
  }

  //! Sets the year.
  /*!
    \param yyyy the year.
  */
  void Date::SetYear(int yyyy)
  {
    year_ = yyyy;
    this->Adjust();
  }

  //! Sets the month.
  /*!
    \param mm the month.
  */
  void Date::SetMonth(int mm)
  {
    month_ = mm;
    this->Adjust();
  }

  //! Sets the day.
  /*!
    \param dd the day.
  */
  void Date::SetDay(int dd)
  {
    day_ = dd;
    this->Adjust();
  }

  //! Sets the hour.
  /*!
    \param hh the hour.
  */
  void Date::SetHour(int hh)
  {
    hour_ = hh;
    this->Adjust();
  }

  //! Sets the minutes.
  /*!
    \param mn minutes.
  */
  void Date::SetMinutes(int mn)
  {
    minutes_ = mn;
    this->Adjust();
  }

  //! Sets the seconds.
  /*!
    \param sc seconds.
  */
  void Date::SetSeconds(double sc)
  {
    seconds_ = sc;
    this->Adjust();
  }

  //! Returns the ordinal number of the day in the year (between 1 and 366).
  /*!
    \return The ordinal number of the day in the year (between 1 and 366).
    \note The ordinal day is often but incorrectly refered to as 'Julian' day.
  */
  int Date::GetOrdinalDay() const
  {
    return this->GetNumberOfDays() + 1;
  }

  //! Returns the number of the day in the year (between 0 and 365).
  /*!
    \return The number of the day in the year (between 0 and 365).
  */
  int Date::GetDayNumber() const
  {
    return this->GetNumberOfDays();
  }

  //! Returns the number of days in the year before the current day.
  /*!
    \return The number of days in the year before the current day. (0 for
    the 1st of January).
  */
  int Date::GetNumberOfDays() const
  {
    int res(0);
    for (int i = 1; i < month_; i++)
      res += month_lengths_[i - 1];
    return res + day_ - 1;
  }

  //! Returns the number of days from a given date.
  /*!
    \param date the reference date.
    \return The number of days between 'date' and the current date
    (positive if the current date is greater than 'date').
  */
  int Date::GetDaysFrom(Date date) const
  {
    int min_year = min(year_, date.GetYear());
    int nb_days(0), nb_days_date(0);
    for (int i = min_year; i < year_; i++)
      nb_days += this->LeapYear(i) ? 366 : 365;
    nb_days += this->GetDayNumber();
    for (int i = min_year; i < date.GetYear(); i++)
      nb_days_date += this->LeapYear(i) ? 366 : 365;
    nb_days_date += date.GetDayNumber();

    return nb_days - nb_days_date;
  }

  //! Returns the number of seconds from a given date.
  /*!
    \param date the reference date.
    \return The number of seconds between 'date' and the current date
    (positive if the current date is greater than 'date').
  */
  double Date::GetSecondsFrom(Date date) const
  {
    int min_year = min(year_, date.GetYear());
    int nb_days(0), nb_days_date(0);
    for (int i = min_year; i < year_; i++)
      nb_days += this->LeapYear(i) ? 366 : 365;
    nb_days += this->GetDayNumber();
    for (int i = min_year; i < date.GetYear(); i++)
      nb_days_date += this->LeapYear(i) ? 366 : 365;
    nb_days_date += date.GetDayNumber();

    return 86400. * double(nb_days - nb_days_date)
      + 3600. * double(hour_ - date.GetHour())
      + 60. * double(minutes_ - date.GetMinutes())
      + double(seconds_ - date.GetSeconds());
  }

  //! Returns the number of hours in the year before the current date.
  /*!
    \return The number of hours in the year before the current date.
  */
  int Date::GetNumberOfHours() const
  {
    return this->GetNumberOfDays() * 24 + hour_;
  }

  //! Returns the number of minutes in the year before the current date.
  /*!
    \return The number of minutes in the year before the current date.
  */
  int Date::GetNumberOfMinutes() const
  {
    return this->GetNumberOfHours() * 60 + minutes_;
  }

  //! Returns the number of seconds in the year before the current date.
  /*!
    \return The number of seconds in the year before the current date.
  */
  double Date::GetNumberOfSeconds() const
  {
    return double(this->GetNumberOfMinutes()) * 60. + seconds_;
  }

  //! Returns the week day.
  /*! Week days are refered as follows: \par
    0 is Monday \par
    1 is Tuesday \par
    2 is Wednesday \par
    3 is Thursday \par
    4 is Friday \par
    5 is Saturday \par
    6 is Sunday \par
    \return The week day number.
  */
  int Date::GetWeekDay() const
  {
    // 1st January 1900 is a Monday.
    // Unknown before.
    if (year_ < 1900)
      return 0;
    int day = 0;
    for (int year = 1900; year < year_; year++)
      day += LeapYear(year) ? 366 : 365;
    day += GetNumberOfDays();
    return day % 7;
  }

  ///////////////
  // OPERATORS //
  ///////////////

  //! Comparison operator <.
  /*!
    \param first_date date.
    \param second_date date.
    \return True if \a first_date is strictly before \a second_date.
  */
  bool operator < (const Date& first_date, const Date& second_date)
  {
    return first_date.GetSecondsFrom(second_date) < 0.;
  }

  //! Comparison operator <=.
  /*!
    \param first_date date.
    \param second_date date.
    \return True if \a first_date is before \a second_date.
  */
  bool operator <= (const Date& first_date, const Date& second_date)
  {
    return first_date.GetSecondsFrom(second_date) <= 0.;
  }

  //! Comparison operator >.
  /*!
    \param first_date date.
    \param second_date date.
    \return True if \a first_date is strictly after \a second_date.
  */
  bool operator > (const Date& first_date, const Date& second_date)
  {
    return first_date.GetSecondsFrom(second_date) > 0.;
  }

  //! Comparison operator >=.
  /*!
    \param first_date date.
    \param second_date date.
    \return True if \a first_date is after \a second_date.
  */
  bool operator >= (const Date& first_date, const Date& second_date)
  {
    return first_date.GetSecondsFrom(second_date) >= 0.;
  }

  //! Comparison operator ==.
  /*!
    \param first_date date.
    \param second_date date.
    \return True if \a first_date is the same date as \a second_date.
  */
  bool operator == (const Date& first_date, const Date& second_date)
  {
    return first_date.GetSecondsFrom(second_date) == 0.;
  }

  //! Comparison operator !=.
  /*!
    \param first_date date.
    \param second_date date.
    \return True if \a first_date is not the same date as \a second_date.
  */
  bool operator != (const Date& first_date, const Date& second_date)
  {
    return first_date.GetSecondsFrom(second_date) != 0.;
  }

  //! Redirection operator <<.
  /*! The date is converted to a string in format "%y-%m-%d %h:%i".
    \param out output stream.
    \param d date to be displayed.
    \return The updated stream.
  */
  ostream& operator << (ostream& out, const Date& d)
  {
    out << d.GetDate("%y-%m-%d %h:%i");
    return out;
  }

}  // namespace Talos.


#define TALOS_FILE_DATE_CXX
#endif
