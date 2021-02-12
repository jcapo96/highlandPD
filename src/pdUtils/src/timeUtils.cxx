#include "timeUtils.hxx"

//*******************************************************
void timeUtils::printTime(const std::string& m){
//*******************************************************  


  gettimeofday(&fTime, NULL);
  double t=fTime.tv_sec+(fTime.tv_usec/1000000.0);

  std::cout << m << ": " << t-fPrevTime << std::endl;
  fPrevTime = t;
}
