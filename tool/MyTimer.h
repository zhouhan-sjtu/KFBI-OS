/*=============================================================================
*   
*   Filename : MyTimer.h
*   Creator : Han Zhou
*   Date : 11/21/21
*   Description : 
*
=============================================================================*/
 
#ifndef _MYTIMER_H
#define _MYTIMER_H


extern "C" {
#include <iostream> 
#include <time.h> 

#ifdef USE_OPENMP
#include <omp.h>
#endif
}

class Timer 
{
#ifdef USE_OPENMP

  public :
    Timer(void) {
      _usedTime = 0.0;
    }

    void start(void) {
      _beginTime = omp_get_wtime(); 
    }

    void stop(void) {
      _endTime = omp_get_wtime(); 
    }

    float getUsedTime(void) {
      _usedTime = _endTime - _beginTime; 
      return _usedTime;
    }

  private :
    float _usedTime; 

    float _endTime; 
    float _beginTime;

#else
  public :
    Timer(void) {
      _usedTime = 0.0;
    }

    void start(void) {
      _beginTime = clock(); 
    }

    void stop(void) {
      _endTime = clock(); 
    }

    float getUsedTime(void) {
      float d = (float)(_endTime - _beginTime); 
      float s = (float)CLOCKS_PER_SEC; 

      _usedTime = d / s; 
      return _usedTime;
    }

  private :
    float _usedTime; 

    time_t _endTime; 
    time_t _beginTime;

#endif

	public :

		void printUsedTime(void) {
			stop();
			float t = getUsedTime();
			std::cout << "CPU time = " << t << " secs." << std::endl;
			start();
		}

}; 

#endif
