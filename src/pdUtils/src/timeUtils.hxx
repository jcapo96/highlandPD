#ifndef timeUtils_h
#define timeUtils_h


#include "pdDataClasses.hxx"
#include <sys/time.h>

class timeUtils{
public:
    
  timeUtils(){fPrevTime=0;}

  void printTime(const std::string& tm);

protected:

  double fPrevTime;
  timeval fTime;
};
 
  
#endif

