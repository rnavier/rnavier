#ifndef RN_TIOBuffer_h
#define RN_TIOBuffer_h

#include <stdarg.h>
#include <iostream>
#include <string>
#include <fstream>
#include <map>

using namespace std ;

//! \brief A Simple class for reading and writing input files and
//! manipulating strings in general
//!
//! This class provides some very
//! simple service routines which are useful for opening, reading and
//! writing input files. 
//!
//! Consider an input file 'foo.txt' in the following form:
//!
//! \verbatim
//! # Comment
//! # Comment
//! [header1] # Comment
//!      5.2 =fValue1 : A given value
//!      #miscellaneous comment
//!      5.3 =fValue2 : A given value
//! [header2]
//!      5 =fValueA : A descriptive string
//!      6 =fValueB : A descriptive string
//! [END]
//! # Comment
//! # Comment
//! [header3]
//!      Name1 =fName1 : A first name
//!      Name2 =fName2 : A second name
//! [END]
//!\endverbatim
//!
//! Calling the following extracts the data from 'foo.txt'
//! \code
//!     TIOBuffer b;
//!     ifstream in("foo.txt") ;
//!
//!     b.read_section(in, "header2") ;
//!     int  vA, vB ;
//!
//!     //This line is error checked the others are not.
//!     if (b.section_has("fValueA"))  {
//!        b.getI("fValueA",vA) ;
//!     }else{
//!        printf("Could not find fValueA!\n") ;
//!     }   
//!     b.getI("fValueB",vB) ;
//!
//!     b.read_section(in, "header1") ;
//!     double v1,v2 ;
//!     b.getD("fValue1", v1) ;
//!     b.getD("fValue2", v2) ;
//!
//!
//!
//! \endcode
//! 
//! \li The sections can be in any order. If the EOF is reached before reaching
//! the section, the file is rewound to the beginning and searced entirely for
//! that section. 
//! \li Once a section is found the section is read until either a
//! new section is found, a [END] is found, or the end of the 
//! file is seen.
//! \li Reading a new section erases the contents of the last section read.  
//! \li All stings are stripped of white space. Thus "header1 " is the same as
//! "header1". The names of headers must be single words 
//! \li Strings which are multiple words separated by white space can not be read
//! with this program.
class TIOBuffer {

   private:
      static const int fNbuffer= 2048 ;
      char fBuffer[fNbuffer] ; 

      void find_name_description(const string &s, string &name, string &description) ;
      std::string fSectionHead ; 
      map <string,string> fSection ;

   public:
      void  readLineI(istream &in, int &i) ;
      void  readLineD(istream &in, double &d) ;
      void  readLineS(istream &in, string &s) ;

      void  skip(istream &in) ;

      void  writeLineI(ostream &out, const int &i, const char *st="") ;
      void  writeLineD(ostream &out, const double &d, const char *st="") ;
      void  writeLineS(ostream &out, const string &s, const char *st="") ;
      void  writeLineB(ostream &out, const bool &b, const char *st="") ;

      char *format(const char *fmt,...) ;

      char *cat(const char *c1, const char *c2) ;

      void read_section(istream &in, const char *head) ;
      void write_section(ostream &out, const char *head) ;
      int  section_has(const char *name) ;
      void getI(const char *name, int &value) ;
      void getD(const char *name, double &value) ;
      void getS(const char *name, string &value) ;
      void getB(const char *name, bool &value) ;
      void strip(string &s) ;


} ;
#endif
