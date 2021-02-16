#include "timeUtils.hxx"


//*******************************************************
timeUtils::timeUtils(size_t ntimes){
//*******************************************************  

  _prevTime=0;
  _times.resize(ntimes,0);  
}

//*******************************************************
void timeUtils::printTime(const std::string& m){
//*******************************************************  


  gettimeofday(&_time, NULL);
  double t=_time.tv_sec+(_time.tv_usec/1000000.0);

  std::cout << m << ": " << t-_prevTime << std::endl;
  _prevTime = t;
}

//*******************************************************
void timeUtils::accumulateTime(size_t index){
//*******************************************************  


  gettimeofday(&_time, NULL);
  double t=_time.tv_sec+(_time.tv_usec/1000000.0);

  _times[index] += t-_prevTime;
  
  _prevTime = t;
}


//*******************************************************
void timeUtils::dumpTimes(){
//*******************************************************  

  std::cout << "times: " << std::endl;
  for (size_t i=0;i<_times.size();i++)
    std::cout << i << " " << _times[i] << std::endl;
  
  
}
