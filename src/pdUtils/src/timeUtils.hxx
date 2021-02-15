#ifndef timeUtils_h
#define timeUtils_h


#include "pdDataClasses.hxx"
#include <sys/time.h>

class timeUtils{
public:
    
  timeUtils(size_t ntimes);

  void printTime(const std::string& tm);

  void accumulateTime(size_t index);

  void dumpTimes();
  
protected:

  double _prevTime;
  timeval _time;

  std::vector<double> _times;

};
 
  
#endif

