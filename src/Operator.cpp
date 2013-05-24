/*
 * Operator.cpp
 *
 *  Created on: 30 Oct 2012
 *      Author: robinsonm
 */

#include "Operator.h"
#include <sstream>
#include <string>
#include <algorithm>
#include "Log.h"

namespace Kairos {
//boost::timer::cpu_timer Operator::global_timer;
double Operator::total_global_time = 0;

Operator::Operator() {
//	timer.stop();
//	global_timer.stop();
	total_time = 0;
	time = 0;
	all_species.clear();

}

int Operator::get_species_index(Species& s) {
	const int n = all_species.size();
	for (int i = 0; i < n; ++i) {
		if (all_species[i] == &s) {
			return i;
		}
	}
	ERROR("did not find species index");
	return -1;
}

void Operator::add_species(Species& s) {
	all_species.push_back(&s);
}

void Operator::resume_timer() {
//	timer.resume();
//	global_timer.resume();
	timer.restart();
}

void Operator::operator ()(const double dt) {
   time += dt;
}

void Operator::reset() {
	time = 0;
}

void Operator::stop_timer() {
	//timer.stop();
	//global_timer.stop();
	const double time = timer.elapsed();
	total_time += time;
	total_global_time += time;
}


template <typename T>
const std::string to_string(const T& data)
{
   std::ostringstream conv;
   conv << data;
   return conv.str();
}

std::string Operator::get_time() {
	//return timer.format();
	return "Time to execute: " + to_string(total_time) + " s (" + get_time_percentage() + ")";
}

std::string Operator::get_global_time() {
	//return global_timer.format();
	return "Time to execute all Operators: " + to_string(total_global_time) + " s";
}

std::string Operator::get_time_percentage() {
//	boost::timer::cpu_times percent;
//	boost::timer::cpu_times this_time = timer.elapsed();
//	boost::timer::cpu_times total = global_timer.elapsed();
//	percent.wall = total.wall / this_time.wall;
//	percent.user = total.user / this_time.user;
//	percent.system = total.system / this_time.system;
//	return boost::timer::format(percent);

	const double percent = 100*total_time/total_global_time;
	return to_string(percent) + "%";
}


}



