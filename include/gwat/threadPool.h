#ifndef THREADPOOL_H
#define THREADPOOL_H
#include <iostream>
#include <functional>
#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

/*! \file 
 *
 * Header file (declarations and definitions because of template functions) for the implementation of a generic thread pool
 */

/*! \brief Default comparator for priority_queue in threadPool -- no comparison
 *
 * First in first out, not sorting
 */
template<class jobtype>
class default_comp
{
public:
	bool operator()(jobtype j, jobtype k)
	{
		return false;
	}
};

/*! \brief Class for creating a pool of threads to asynchronously distribute work
 *
 * Template parameters: 
 *
 * jobtype defines a structure or class that represents a job or task 
 *
 * comparator defines how to compare jobs for sorting the list
 *
 * Default options correspond to jobs being defined by an integer job_id, and no sorting of the list (first in first out)
 */
template<class jobtype=int, class comparator=default_comp<jobtype>>
class threadPool
{
public: 

/*! \brief Constructor -- starts thread pool running
 */
	explicit threadPool(std::size_t numThreads, std::function<void(int, jobtype)> work_fn)
	{
		work_fn_internal = work_fn;
		numThreads_internal = numThreads;

		start(numThreads);
	}
/*! \brief Destructor -- stops threads
 */
	~threadPool()
	{
		stop();	
	}
	
/*! \brief Places jobs in queue to wait for scheduling
 *
 * job_id is sorted if a comparator is provided
 */
	void enqueue(jobtype job_id)
	{
		{
			std::unique_lock<std::mutex> lock{EventMutex};
			Tasks.emplace(std::move(job_id));
		}
		EventVar.notify_one();

	}
/*!\brief Get the number of threads being used by the thread pool
 */
	int get_num_threads()
	{
		return numThreads_internal;
	}
/*! \brief Get the current length of the job queue
 */
	int get_queue_length()
	{
		return Tasks.size();
	}
	
	
private:
/*! Vector for thread ids*/
	std::vector<std::thread> Threads;

/*! Condition Variable*/
	std::condition_variable EventVar;

/*! Lock to prevent memory races*/
	std::mutex EventMutex;

/*! Boolean stop condition*/
	bool stopping=false;

/*! Number of threads in pool*/
	int numThreads_internal;

/*! Function for each thread to perform -- Takes arguments thread_id and jobtype job*/
	std::function<void(int, jobtype)> work_fn_internal;

/*! Queue of jobs*/
	std::priority_queue<jobtype,std::vector<jobtype>, comparator> Tasks;
	
/*! \brief Starts thread pool -- launches each thread, which continually check for work
 */
	void start(std::size_t numThreads)
	{
		for(auto i =0u; i<numThreads; i++)
		{
			Threads.emplace_back([=, this]{
				while(true)
				{
					jobtype j;
					{
						std::unique_lock<std::mutex> lock{EventMutex};

						EventVar.wait(lock,[=, this]{return stopping || !Tasks.empty(); });
						
						if (stopping && Tasks.empty())
							break;	
						j = std::move(Tasks.top());
						Tasks.pop();
					}
					work_fn_internal(i, j);
					
				}
			});
		}
	}
/*! \brief Stops thread pool 
 *
 * Waits for all threads to end and joins them
 *
 * Finishes the thread pool before ending
 */
	void stop() noexcept
	{
		//std::cout<<std::endl;
		//std::cout<<"Stop initiated -- waiting for threads to finish"<<std::endl;
		{
			std::unique_lock<std::mutex> lock{EventMutex};
			stopping = true;
		}
		
		EventVar.notify_all();
		
		for(auto &thread: Threads)
			thread.join();
	}
};

#endif
