#include <iostream>
#include <cstdlib>
#include "RNavierPath.h"

using namespace std;

static string RNavierNametag = "";
static string RNAVIERDATA= ""  ;

const string &GetRNAVIERDATA() 
{
   static bool first = true;
   // Get the path from the environment variable RNAVIERDATA.  Append a "/" if
   // not already there.  Assign the gRNAVIERDATA to the result
   if (first) {
      char *path =getenv("RNAVIERDATA")  ;
      if (!path)  {
         cerr << "** GetRNAVIERDATA ** Warning environment variable RNAVIERDATA not set" << endl;
         RNAVIERDATA = "" ;
      } else {
         string s = path;
         if (s.rbegin()!=s.rend()) {
            char c = *s.rbegin() ;
            if (c!='/') s+= "/" ;
         } 
         RNAVIERDATA = s;
      }
      first = false ;
   }
   return RNAVIERDATA ;
}
const string &GetRNavierNametag() 
{
   if (RNavierNametag.size() == 0)  {
      cerr << "*** GetRNavierNametag() *** The current value of nametag is not set. Call SetRNavierNametag() before calling this function" << endl;
      exit(EXIT_FAILURE) ;
   }
   return RNavierNametag ;
}
void SetRNavierNametag(const string &nametag) 
{
   if (nametag.size() == 0)  {
      cerr << "*** gSetNametag() *** Bad input -- nametag is empty! " << endl;
      exit(EXIT_FAILURE) ;
   }
   RNavierNametag  = nametag;
}



