/* Library */
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <thread>
#include <time.h>
#include <unordered_map>
#include <string>
#include <numeric>
#include <algorithm>
#include <iomanip>


/* My library */
#include "cellular-automata.hpp"
#include "cash-display.hpp"

/* Other headers */
#include "automaton.hpp"

#ifndef RUNTIME
#define RUNTIME

std::vector<unsigned> runtime(std::unordered_map<std::string, unsigned>& s_ops, std::unordered_map<std::string, double>& a_prms, std::string experiment_tag, std::vector<unsigned>& populations, std::vector<std::vector<Automaton>>& virions, std::vector<double> influx_populations);

unsigned draw_gap(unsigned time, unsigned pixel_visit, unsigned dt, unsigned pop_count); 

#endif 