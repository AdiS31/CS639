#pragma once

#include <chrono>
#include <cstring>
#include <iostream>

struct Timer
{
    using clock_t = std::chrono::high_resolution_clock;
    using time_point_t = std::chrono::time_point<clock_t>;

    time_point_t mStartTime;
    time_point_t mStopTime;
	double runtime;
	
    void Start()
    {
        mStartTime = clock_t::now();
		runtime = 0.0;
    }
	
	void Stop()
	{
		mStopTime = clock_t::now();
		std::chrono::duration<double, std::milli> elapsedTime = mStopTime - mStartTime;
		runtime = elapsedTime.count();
	}
	
    void Stop(const std::string& msg)
    {
        mStopTime = clock_t::now();
        std::chrono::duration<double, std::milli> elapsedTime = mStopTime - mStartTime;
		runtime = elapsedTime.count();
        std::cout << "[" << msg << elapsedTime.count() << "ms]" << std::endl;
    }
	
	double getRuntime()
	{
		return runtime;
	}
    
};
