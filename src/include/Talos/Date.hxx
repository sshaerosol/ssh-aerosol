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


#ifndef TALOS_FILE_DATE_HXX


#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>


namespace Talos
{

  using namespace std;

  //! Class Date.
  class Date
  {
  private:
    //! Year.
    int year_;
    //! Month.
    int month_;
    //! Day.
    int day_;
    //! Hour.
    int hour_;
    //! Minutes.
    int minutes_;
    //! Seconds.
    double seconds_;
    //! Month lengths;
    vector<int> month_lengths_;

    void LeapYearAdjust();
    bool IsValid();
    void Adjust();

  public:
    Date();
    Date(const Date& date);
    Date(string date);
    Date(int yyyymmdd);
    Date(int yyyy, int mm, int dd = 1,
         int hh = 0, int mn = 0, double sc = 0);

#ifndef SWIG
    Date& operator=(const Date&);
    Date& operator=(string date);
#endif
    void SetDate(string date);
    void SetDate(int yyyymmdd);
    void SetDate(int yyyy, int mm, int dd = 1,
                 int hh = 0, int mn = 0, double sc = 0);

    bool LeapYear(int year) const;
    bool LeapYear() const;

    int GetDate() const;

    string GetDate(const string& format) const;
    template <class T>
    void GetDate(const string& format, T& date) const;

    int GetYear() const;
    int GetMonth() const;
    int GetDay() const;
    int GetHour() const;
    int GetMinutes() const;
    double GetSeconds() const;

    void AddYears(int nb_yy);
    void AddMonths(int nb_mm);
    void AddDays(int nb_dd);
    void AddHours(int nb_hh);
    void AddMinutes(int nb_mn);
    void AddSeconds(double nb_sc);

    void SetYear(int yyyy);
    void SetMonth(int mm);
    void SetDay(int dd);
    void SetHour(int hh);
    void SetMinutes(int mn);
    void SetSeconds(double sc);

    int GetOrdinalDay() const;
    int GetDayNumber() const;
    int GetNumberOfDays() const;
    int GetDaysFrom(Date date) const;
    double GetSecondsFrom(Date date) const;

    int GetNumberOfHours() const;
    int GetNumberOfMinutes() const;
    double GetNumberOfSeconds() const;

    int GetWeekDay() const;
  };

  // Comparisons.
  bool operator < (const Date& first_date, const Date& second_date);
  bool operator <= (const Date& first_date, const Date& second_date);
  bool operator > (const Date& first_date, const Date& second_date);
  bool operator >= (const Date& first_date, const Date& second_date);
  bool operator == (const Date& first_date, const Date& second_date);
  bool operator != (const Date& first_date, const Date& second_date);

#ifndef SWIG
  // Redirection.
  ostream& operator << (ostream& out, const Date& d);
#endif

}  // namespace Talos.


#define TALOS_FILE_DATE_HXX
#endif
