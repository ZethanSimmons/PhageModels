/* Library */
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <thread>
#include <time.h>
#include <string>
#include <unordered_map>
#include <math.h>

/* My library */
#include "cellular-automata.hpp"
#include "cash-display.hpp"

/* Other headers */
#include "automaton.hpp"
#include "runtime.hpp"

int main(int argc, char* argv[])
{   
    //initialize arguments to be fed to runtime later, from slurm script

    //population is used to initialize the dominant population in "initial_populations" later
    int dom_population = 000; 
    int perturbation_type = 0x000;
    std::string output_location = ".";
    std::unordered_map<std::string, unsigned> s_ops = {}; //"simulation options" - integer-basedd options to feed into the simulation. Some of them are actually parameters for the automatons
    std::unordered_map<std::string, double> a_prms = {}; //"automaton parameters" - decimals feed into runtime as simulation parameters, mostly applicable to automatons
    std::string experiment_tag = "_ex1";

    //SIMULATION RUN OPTIONS
    s_ops["diffuse_switch"] = false; //not implemented yet
    s_ops["random_mixing"] = false; //so far, random mixing has been hard-coded into runtime.cpp
    s_ops["grid_dimension"] = 200; // used for width and height of the grid, in pixels/dots
    s_ops["screenheight"] = 2000; //has no impact on simulation - used only for live simulation display
    s_ops["the_seed"] = 1234567; //used to seed rng
    s_ops["max_time"] = 10000; //maximum time to run the simulation
    //end_early: 0 = false; 1 = end early if there is only one type of bacteria left, considering plasmids as well; 2 - end early if there is only one type of chromosomes left, regardless of plasmids 
    s_ops["end_early"] = 0; //boolean to end simulation early if certain conditions are met (see runtime function)

    s_ops["perturbation_magnitude"] = 0; //controls the size of a single population perturbation put in place at the start of the simulation
    s_ops["evolve_parameters_switch"] = 0; //controls whether to evolve parameters of the simulation

    unsigned* grid_dimension = &s_ops["grid_dimension"]; //can be used in calculations for initial condition setup in main(), together with gridsize() lambda
    auto gridsize = [&grid_dimension]() { return *grid_dimension * *grid_dimension; };

    //OUTPUT OPTIONS
    s_ops["output_data_flag"] = true; //outputs population data at given intervals of time -- see output_gap
    s_ops["output_gap"] = 200; //time interval to output population data

    s_ops["show_display_flag"] = false; //boolean to show live image output
    s_ops["draw_gap"] = gridsize(); //within a unit of time, interval to make live image output
    s_ops["draw_gap_time"] = 100; //interval for decide which timesteps to make live image output

    s_ops["make_movie_flag"] = false; //outputs image of grid as jpgs, at chosen time interval
    s_ops["micro_movie_gap"] = 100; //outputs movies close together at the very start of the simulation
    s_ops["micro_movie_limit"] = 1000;
    s_ops["init_movie_gap"] = 1000;
    s_ops["init_movie_limit"] = 10000;
    s_ops["movie_gap"] = 100000; //time interval to make jpg image outputs

    s_ops["test_flag"] = false; //can activate certain testing behaviours of runtime() function
    s_ops["delay_period"] = 0; //if test_flag is true, can delay each live image output for a chosen amount of ms



    //AUTOMATON PARAMETERS

    a_prms["delta_t"] = 1.0; //can be used to slow model down, if set less than 1.0
    a_prms["v_diffusion"] = 0.06; //diffusion rate of virions
    a_prms["b_growth"] = 0.3; //background growth rate of bacteria

    // bacterial death rate, universal antibiotics, and l rate must sum to 1 or less
    a_prms["b_death"] = 0.06; //background deathrate of bacteria
    a_prms["conjugation"] = 0.14; //conjugation rate of plasmids
    a_prms["conjugation_010"] = 0.14; //conjugation rate of plasmids when bacteria is infected by 010 plasmids
    a_prms["s_loss"] = 0.0015; //segregative loss rate of plasmids during bacterial reproduction

    a_prms["b_resistance"] = 0.3; //strength of antibiotic resistance carried by bacteria
    a_prms["p_resistance"] = 0.3; //strength of antibiotic resistance carried by plasmids
    a_prms["v_resistance"] = 0.3; //strength of antibiotic resistance carried by phages

    a_prms["p_cost"] = 0.00225; //growth cost of plasmids - subtracted from b_growth
    a_prms["r_cost"] = 0.0015; //growth cost of resistance - subtracted once from b_growth for each carrier of antibiotic resistance present in a given bacteria (chromosome, plasmid, phage are possible)

    //PHAGE TRAITS
    s_ops["burst_size"] = 4; //number of virions produced upon viral lysis of bacteria
    s_ops["burst_size_001"] = 4; //burst_size of bacterial infexted by 001 phages
    a_prms["lysogenization_rate"] = 0.0; //rate of lysogenization after infection of bacteria by phages - non-lysogenized, infected bacteria are lysed
    a_prms["lysogenization_rate_001"] = 0.0; //rate of lysogenization after infection of bacteria by 001 phages
    a_prms["lysis_rate"] = 0.00; //background lysis rate of lysogens
    a_prms["lysis_rate_001"] = 0.00; //background lysis rate of bacteria
    a_prms["v_degradation_rate"] = 0.06; //background degredation rate of virions
    a_prms["infection_rate"] = 0.06; //background rate of infection by a single virion into a bacteria at the same location on the grid
    a_prms["v_loss_rate"] = 0.00015; //loss rate of phages from lysogens - for simplicity, treated as a segregative loss rate, like phages, and occurs at reproduction, even though in real life viral loss can occur outside of reproduction 
    

    s_ops["death_lysis_switch"] = 0; //controls whether lysogens are lysed when dying
    
    //ANTIBIOTIC CONFIGURATION
    // set antibiotic_time_periodicity and antibiotic_pixel_periodicity less than the block parameters, respectively, to ensure constant antibiotics
    
    a_prms["antibiotics"] = 0.3; //antibiotic-associated deathrate of bacteria - additional to b_death
    s_ops["antibiotic_time_periodicity"] = 100; //time interval after which environmental antibiotics is switched on
    s_ops["antibiotic_time_limit"] = 1; //time interval after which environmental antibiotics is switched off
    s_ops["antibiotic_pixel_periodicity"] = 1; //pixel-visit interval after which environmental antibiotics is switched on. Within a single timestep, there are gridsize() pixelvisits
    s_ops["antibiotic_pixel_limit"] = 2; //within a single timestep, pixel-visit interval after which environmental antibiotics is switched off. There are gridsize() pixelvisits per timestep
    

    std::vector<double> experimental_influx = { //used later during main simulation to control influx of bacteria
    0.0, 
    0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, //no phage
    0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, //susceptible phage
    0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000  //AB-resistant phage
    }; 

    std::unordered_map<std::string, int> bacterial_types = {
    {"000", 0},  
    {"100", 1}, {"110", 2}, {"190", 3}, {"900", 4}, {"910", 5}, {"990", 6}, //no phage
    {"101", 7}, {"111", 8}, {"191", 9}, {"901", 10}, {"911", 11}, {"991", 12}, //susceptible phage
    {"109", 13}, {"119", 14}, {"199", 15}, {"909", 16}, {"919", 17}, {"999", 18}  //AB-resistant phage
    }; 

    //unpack argv arguments into the 5 arguments for runtime()
    std::string key;
    for (int i = 1; i < argc; i += 2) {
        if (std::string(argv[i]).substr(0,4) == "--a_") {
            key = std::string(argv[i]).substr(4, std::string(argv[i]).size() - 4);
            std::cout<<key<<"\n";
            a_prms[key] = std::stod(argv[i+1]);
        } else if (std::string(argv[i]).substr(0,4) == "--s_") {
            key = std::string(argv[i]).substr(4, std::string(argv[i]).size() - 4);
            std::cout<<key<<"\n";
            s_ops[key] = std::stoi(argv[i+1]);
        } else if (std::string(argv[i]).substr(0,4) == "--i_") {
            key = std::string(argv[i]).substr(4, std::string(argv[i]).size() - 4);
            std::cout<<key<<"\n";
            experimental_influx[bacterial_types[key]] = std::stod(argv[i+1]);
        } else if (std::string(argv[i]) == "--population") {
            dom_population = std::stoi(argv[i+1], nullptr, 16);
        } else if (std::string(argv[i]) == "--perturbation_type") {
            perturbation_type = std::stoi(argv[i+1], nullptr, 16);
        } else if (std::string(argv[i]) == "--output_location") {
            output_location = std::string(argv[i+1]);
        } else {
            std::cerr << "Error: Unrecognized argument" << argv[i] << std::endl;
        }
    }

    std::cout << "Finished reading arguments\n";

    //initialize different populations of bacteria, immigrating bacteria, and virions on grids, to feed into the runtime() function
    std::vector<unsigned> initial_populations(gridsize()); 
    std::vector<double> no_influx = {
        0.0,                          //empty
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //no phage
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //susceptible phage
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0  //AB-resistant phage
        }; 
  
    std::vector<std::vector<Automaton>> empty_virions(gridsize(), std::vector<Automaton>(0)); 

    std::cout << "Finished null population setup\n";

    //setup rng
    std::seed_seq seed({s_ops["the_seed"]});
    std::mt19937_64 random(seed);
    std::uniform_int_distribution<int> rand_grid(0,gridsize());

    std::cout << "Finished rng setup\n";

    //setup output locations to output movies of equilibriation run and true experimental run
    std::string loc1 = output_location + "/00_experiment_run_";
    std::string loc2 = output_location + "/01_experiment_run_evolved_";

    //pop_labels is later used to construct output1 and output2 outuput filenames
    std::vector<std::string> pop_labels = {
        "empty", 
        "c00", "cp0", "cP0", "C00", "Cp0", "CP0", //no phage
        "c0v", "cpv", "cPv", "C0v", "Cpv", "CPv", //susceptible phage
        "c0V", "cpV", "cPV", "C0V", "CpV", "CPV"  //AB-resistant phage
        }; 


    //setup 
    auto dom_pop_it = std::max_element(initial_populations.begin(), initial_populations.end());
    auto dom_pop_index = std::distance(initial_populations.begin(), dom_pop_it);
    fill(initial_populations.begin(), initial_populations.end(), dom_population);

    //Step 1: Run the initial 10,000-timestep run to determine the time 0 starting layout
    std::string data_info = "";
    std::string output1 = loc1 + "_CA";
    std::string output2 = loc2 + "_CA";


    if (s_ops["run1"])
    {
    // Step 1: Setup first 10,000,000-timestep run with no influx, to evolve parameters
    s_ops["max_time"] = 10000001; // Set max_time for the 1,000,000-timestep run
    s_ops["evolve_parameters_switch"] = true; // Enable parameter evolution
    s_ops["output_gap"] = 200;

    // Step 2: Run the initial 10,000,000-timestep run to determine the time 0 starting layout
    std::cout << "Start initial 10,000,000-timestep run\n";
    std::vector<unsigned> initial_results = runtime(s_ops, a_prms, output1, initial_populations, empty_virions, no_influx);
    std::cout << "Finished initial 10,000,000-timestep run\n";
    }

    if (s_ops["run2"])
    {
    // Step 6: Run the 3,000,000-timestep run with evolve_parameters_switch set to true and the evolved rates
    s_ops["evolve_parameters_switch"] = true; // Enable parameter evolution
    s_ops["max_time"] = 10000001; // Set max_time for the 1,000,000-timestep run
    s_ops["output_gap"] = 200;

    std::cout << "Start 10,000,000-timestep run with evolve_parameters_switch = true\n";
    runtime(s_ops, a_prms, output2, initial_populations, empty_virions, experimental_influx);
    std::cout << "Finished 10,000,000-timestep run\n";
    }
}

