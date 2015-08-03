#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <cctype>
#include <sstream>
#include <string.h>
#include "TIOBuffer.h"

using namespace std ;

//! Read a c++ style string from the input output buffer
void TIOBuffer::readLineS(istream &in, string &s) 
{
   skip(in) ;
   in >> s ;
   in.getline(fBuffer,fNbuffer) ;
}


//! Read an integer int &i from an input file.  See readLineD.
void TIOBuffer::readLineI(istream &in, int &i) 
{
   skip(in) ;
   in >> i ;
   in.getline(fBuffer,fNbuffer) ;
}

//! Read in a double from an  input file.

//! Skip over all lines beggining with # and white space,
//! read in a double, and then skip over trailing characters 
//! in that line.
//! For example, an input file that looks like  this
//! \verbatim
//! #line 1 
//! #line 2 
//!      5.2  a variable
//!\endverbatim
//! Calling  
//! \code 
//! b.readLineD(in, xmin) ; 
//! \endcode 
//! would read in 5.2 into xmin, skipping over the first and second
//! line and the variable description.
void TIOBuffer::readLineD(istream &in, double &d) 
{
   
   skip(in) ;
   in >> d ;
   in.getline(fBuffer,fNbuffer) ;

}

//! Returns a string formatted with c-style formats.

//! This routine is useful for printing columns.    For example,
//! calling:  
//! \code clog << b.format("%10.5 %10.5f\n", 5.2, 5.2) << endl; \endcode 
//! produces:
//!\verbatim   5.20000    5.20000 \endverbatim
char *TIOBuffer::format(const char *fmt,...)
{
   va_list ap ;
   va_start(ap, fmt) ;
   int n = vsnprintf(fBuffer, fNbuffer, fmt, ap);
   // old vsnprintf's return -1 if string is truncated new ones return
   // total number of characters that would have been written
   if (n == -1 || n >= fNbuffer) {
      cerr << "** TIOBuffer::format ** Exceeded buffer length " <<endl;
   }

   va_end(ap) ;
   return fBuffer ;
}

//! Skips over white space and all lines beginning with '#' until the next input is seen.
void TIOBuffer::skip(istream &in) 
{
   char c ;
   
   while(1) {
      in.get(c) ;
      if (!isspace(c)) {
         if (c=='#') {
            in.getline(fBuffer,fNbuffer) ;
         } else if (c=='[') {
            in.getline(fBuffer,fNbuffer) ;
         } else{
            in.putback(c) ;
            break;
         }
      }
   }
}

//! \brief Write an integer to an output file together with descriptive string 

//! Example: 
//!\code 
//!     b.writeLine(clog, 400, "fImax ; A description"); 
//!\endcode returns: 
//!\verbatim
//!        fImax = 400 ; A description \endverbatim
void  TIOBuffer::writeLineI(ostream &out, const int &i, const char *st) 
{
   string name, description ;
   find_name_description(st, name, description) ;
   out << format("%-20s = %-12d ; %s\n", name.c_str(), i, description.c_str()) ;  
}

//! Writes a double to an output file, see  writeLineI.
void  TIOBuffer::writeLineD(ostream &out, const double &d, const char *st) 
{
   string name, description ;
   find_name_description(st, name, description) ;
   out << format("%-20s = %-12g ; %s\n", name.c_str(), d, description.c_str()) ;  
}
//! Writes a string to an output file, see  writeLineI
void  TIOBuffer::writeLineS(ostream &out, const string &s, const char *st) 
{
   string name, description ;
   find_name_description(st, name, description) ;
   out << format("%-20s = %-12s ; %s\n", name.c_str(), s.c_str(), description.c_str()) ;  
}
//! Writes a string to an output file, see  writeLineI
void  TIOBuffer::writeLineB(ostream &out, const bool &b, const char *st) 
{
   string name, description ;
   find_name_description(st, name, description) ;
   if  (b) {
     out << format("%-20s = %-12s ; %s\n", name.c_str(), "True", description.c_str()) ;  
   } else {
     out << format("%-20s = %-12s ; %s\n", name.c_str(), "False", description.c_str()) ;  
   }
}

//! Returns a pointer to a sctring containing the catenation of string c1 and string c2.
char *TIOBuffer::cat(const char *c1, const char *c2) 
{
   strcpy(fBuffer, c1) ; 
   strcat(fBuffer, c2) ;
   return fBuffer ;
}


//! Read an input file section.
void TIOBuffer::read_section(istream &in, const char *head) 
{

   string heads(head) ;
   fSectionHead = heads ;

   string s1 ;

   int is_first_pass = 1 ;
   char c ;
   fSection.clear() ;
   while(1)  {
      //v = in.peek() ;
      in.peek() ;
      if ( (in.eof()) && (is_first_pass) ) {
         is_first_pass = 0 ;
         in.clear() ;
         in.seekg(0, ios::beg)  ; 
         continue ;
      } else if ( (in.eof()) && (!is_first_pass) ) {
         clog <<"** TIOBuffer::read_section ** Section " << heads 
              << " is not found" << endl;
         exit(-1) ;
         return ;
      }
      in.get(c) ;
      if (!isspace(c)) {
         if (c=='#') {
            in.getline(fBuffer,fNbuffer) ;
         } else if (c=='[') {
            in.get(fBuffer, fNbuffer, ']') ;
            s1 = fBuffer ;
            strip(s1) ;
            in.getline(fBuffer, fNbuffer) ;
            if (s1 == heads) break ;
         } else {
            in.getline(fBuffer, fNbuffer) ;
         }
      }
   }

   while(1) {
      in.peek() ;
      if (in.eof())  {
         in.clear() ;
         in.seekg(0, ios::beg)  ; 
         break ;
      }

      in.get(c)  ;
      if (!isspace(c)) {
         if (c == '#') {
            in.getline(fBuffer, fNbuffer) ;
            continue ;
         }else if (c == '[') {
            //This is the start of a new fSection
            in.putback(c) ;
            break ;
         }else {
            in.putback(c) ;
         }

         in.getline(fBuffer, fNbuffer) ;
         string line = fBuffer ;
         string::size_type i1,i2,i3 ; 

         // Strip off trailing comments
         i1 = line.find(";")  ;
         i2 = line.find("#")  ;
         i3 = (i1>= i2 ? i1 : i2) ;
         string line2 = line.substr(0,i3) ;

         if ((i1=line2.find("="))!=string::npos) {
            line2.replace(i1,1," ") ;
         } else {
            cerr << "** TIOBuffer::read_section ** Error in reading line = \n" << line << endl;
            cerr << "Line must have the form " <<endl;
            cerr << "Variable  = value ; Description" << endl;
            cerr << "Did you forget the equal sign? " << endl;
            exit(-1) ;
         }
         istringstream ist(line2) ;
         string value, name  ;
         ist >> name  >> value ;
         fSection.insert(make_pair(name,value)) ;
      }
   }
   return ;
}

//! Returns 1/0 if the last read section has a member named, name.
int  TIOBuffer::section_has(const char *name) 
{
   string s = name ;
   map<string,string>::iterator p  = fSection.find(s) ;

   if (p==fSection.end()) {
      return 0 ;
   }else{
      return 1 ;
   }
}
//! Gets an integer named name from the last read section.
void TIOBuffer::getI(const char *name, int &value) 
{
   string s = name ;

   map<string,string>::iterator p  = fSection.find(s) ;
   if (p == fSection.end()) { 
      cerr << "** TIOBuffer::getI ** name -- " << name <<
         " -- was not found in section " << fSectionHead << endl;
      exit(EXIT_FAILURE) ;
      return ;
   }

   string v = fSection[name] ;
   value = atoi(v.c_str()) ;
}
//! Gets a double from the last read section.
void TIOBuffer::getD(const char *name, double &value) 
{
   string s = name ;
   map<string,string>::iterator p  = fSection.find(s) ;

   if (p == fSection.end()) { 
      cerr << "** TIOBuffer::getD ** name -- " << name <<
         " -- was not found in section " << fSectionHead << endl;
      exit(EXIT_FAILURE) ;
      return ;
   }

   string v = fSection[name] ;
   value =  atof(v.c_str()) ;
}
//! Gets a boolean from the last read section
void TIOBuffer::getB(const char *name, bool &value) 
{
   string s = name ;
   map<string,string>::iterator p  = fSection.find(s) ;

   if (p == fSection.end()) { 
      cerr << "** TIOBuffer::getB ** name -- " << name <<
         " -- was not found in section " << fSectionHead << endl;
      exit(EXIT_FAILURE) ;
      return ;
   }
   string v = fSection[name] ;
   if (v=="True"||v=="true"||v=="on"||v=="ON"||v=="On"||v=="Yes"||v=="yes"||v=="1") {
      value = true ;
   } else {
      value = false ;
   }
}
//! Gets a string from the last read section.
void TIOBuffer::getS(const char *name, string &value) 
{
   string s = name ;
   map<string,string>::iterator p  = fSection.find(s) ;

   if (p == fSection.end()) { 
      cerr << "** TIOBuffer::getS ** name -- " << name <<
         " -- was not found in section " << fSectionHead << endl;
      exit(EXIT_FAILURE) ;
      return ;
   }

   value = fSection[name] ;
}


//! Strips leading and trailing spaces from a string.
void TIOBuffer::strip(string &s) 
{
   unsigned int i1, i2, i3, i4 ;
   unsigned int i ;

   i1 = 0 ;
   i2 = 0 ;
   for (i = 0 ; i < s.size() ; i++)  {
      if (!isspace(s[i])) {
         i2 = i ;
         break ;
      }
   }
   if (i2 !=0 ) {
      s.erase(0,i2) ;
   }
   
   i3 = s.size() ;
   i4 = s.size() ;
   for (i = i2 ; i < s.size() ; i++) {
      if (isspace(s[i])) {
         i3 = i ; 
         break ;
      }
   }
   if (i3 != s.size() ) {
      s.erase(i3,s.size()) ;
   }
}

void TIOBuffer::write_section(ostream &out, const char *head) 
{
   out << "[" << head <<"]\n" ;
   map<string,string>::iterator p  ;

   for (p=fSection.begin() ; p != fSection.end(); ++p) {
      out << p->first <<" = "<< p->second << "\n" ;
   }
}
void TIOBuffer::find_name_description(const string &s, string &name, string &description) 
{
   string line = s ;
   string::size_type i1,i2,i3 ;

   // Strip off trailing comments
   i1 = line.find(";")  ;
   i2 = line.find("#")  ;
   i2 = (i1<= i2 ? i1 : i2) ;
   i3 = line.find(":")  ;
   i3 = (i2<= i3 ? i2 : i3) ;

   name = line.substr(0,i3) ;
   strip(name) ;

   if (i3!=string::npos) {
      description = line.substr(i3+1,string::npos) ;
   }
}


   
