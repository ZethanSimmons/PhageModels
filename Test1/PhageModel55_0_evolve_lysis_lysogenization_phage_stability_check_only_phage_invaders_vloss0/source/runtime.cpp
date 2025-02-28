#include "runtime.hpp"

unsigned testing(std::string input) {
  // std::cout << input << "\n";
  return 0;
};

unsigned testing(int input) {
  // if ((input & 0x005 == 0x005) || (input & 0x00D == 0x00D))
  //  std::cout << input << "\n";
  return 0;
};

void make_output(){

}

void initialize_grid() {

}

std::vector<unsigned> runtime(std::unordered_map<std::string, unsigned>& s_ops, std::unordered_map<std::string, double>& a_prms, std::string experiment_tag, std::vector<unsigned>& populations, std::vector<std::vector<Automaton>>& virions, std::vector<double> influx_populations)
{ 
    //BLOCK 1 - UNPACK ARGUMENTS

    //s_ops - simulation run options; options that are integers
    unsigned test_flag = s_ops["test_flag"];
    unsigned output_data_flag = s_ops["output_data_flag"];
    unsigned delay_period = s_ops["delay_period"];
    unsigned rows_n = s_ops["grid_dimension"];
    unsigned cols_n = s_ops["grid_dimension"];
    unsigned gridsize = rows_n * cols_n;
    unsigned max_time = s_ops["max_time"];
    unsigned scale = s_ops["screenheight"] / rows_n / 2;
    unsigned draw_gap = s_ops["draw_gap"];
    unsigned movie_gap = s_ops["movie_gap"]; //gap for image output throughout simulation
    unsigned init_movie_gap = s_ops["init_movie_gap"]; //gap for image output up to movie_lim
    unsigned init_movie_lim = s_ops["init_movie_limit"]; //limit for init_movie_gap to apply for image output
    unsigned micro_movie_gap = s_ops["micro_movie_gap"]; //gap for image output up to micro_movie_lim
    unsigned micro_movie_lim = s_ops["micro_movie_limit"]; //limit for micro_movie_gap to apply for image output
    unsigned output_gap = s_ops["output_gap"];
    unsigned the_seed = s_ops["the_seed"];
    unsigned draw_gap_time = s_ops["draw_gap_time"];
    unsigned end_early = s_ops["end_early"];

    // set antibiotic_time_periodicity and antibiotic_pixel_periodicity less than the block parameters, respectively, to ensure constant antibiotics; 
    unsigned antibiotic_time_periodicity = s_ops["antibiotic_time_periodicity"];
    unsigned antibiotic_time_limit = s_ops["antibiotic_time_limit"];
    unsigned antibiotic_pixel_periodicity = s_ops["antibiotic_pixel_periodicity"];
    unsigned antibiotic_pixel_limit = s_ops["antibiotic_pixel_limit"];

    //a_prms - automaton parameters; mostly model parameters, which are doubles 
    Automaton::set_delta_t(a_prms["delta_t"]);
    double delta_t = a_prms["delta_t"];
    double v_diffusion = a_prms["v_diffusion"] * delta_t;
    double b_growth = a_prms["b_growth"] * delta_t;
    double b_death = a_prms["b_death"] * delta_t; // bacterial death rate must sum to 1 or less with plasmid conjugation_rate and antibiotics
    double conjugation = a_prms["conjugation"];
    double conjugation_010 = a_prms["conjugation_010"];
    double conjugation_t0 = a_prms["conjugation_t0"];
    double s_loss = a_prms["s_loss"]; // segregation loss
    double influx_rate = std::accumulate(influx_populations.begin(),influx_populations.end(),0.0);

    double b_resistance = a_prms["b_resistance"];
    double p_resistance = a_prms["p_resistance"];
    double v_resistance = a_prms["v_resistance"];

    double volatile influx_chance;

    double p_cost = a_prms["p_cost"];
    double r_cost = a_prms["r_cost"];

    //virion configuration
    unsigned burst_size = s_ops["burst_size"];
    unsigned burst_size_001 = s_ops["burst_size_001"];
    double lysogenization_rate = a_prms["lysogenization_rate"];
    double lysogenization_rate_001 = a_prms["lysogenization_rate_001"];
    double lysogenization_rate_t0 = a_prms["lysogenization_rate_t0"];
    double lysis_rate = a_prms["lysis_rate"];
    double lysis_rate_001 = a_prms["lysis_rate_001"];
    double lysis_rate_t0 = a_prms["lysis_rate_t0"];
    double v_degradation_rate = a_prms["v_degradation_rate"];
    double infection_rate = a_prms["infection_rate"];
    double v_loss_rate = a_prms["v_loss_rate"];

    

    double antibiotic_pressure = a_prms["antibiotics"] * delta_t;
    bool antibiotic_trigger;

    bool diffuse_switch = s_ops["diffuse_switch"]; 
    bool death_lysis_switch = s_ops["death_lysis_switch"];
    bool random_mixing = s_ops["random_mixing"];
    bool show_display_flag = s_ops["show_display_flag"];
    bool make_movie_flag = s_ops["make_movie_flag"];
    bool evolve_parameters_switch = s_ops["evolve_parameters_switch"];

    //BLOCK 2 - SETUP SIMULATION GRID

    CA2D<Automaton>* ca_grid = new CA2D<Automaton>(rows_n, cols_n);
    CA2D<std::vector<Automaton>>* virion_grid = new CA2D<std::vector<Automaton>>(rows_n,cols_n);


    //setup reference automata representing different types of bacteria, etc, to use to initialize CA
    std::vector<unsigned> bacteria_types = {
      0, 
      0x100, 0x110, 0x190, 0x900, 0x910, 0x990,
      0x101, 0x111, 0x191, 0x901, 0x911, 0x991,
      0x109, 0x119, 0x199, 0x909, 0x919, 0x999
    };
    //reference populations is used to initialize the populations on the grid
    std::unordered_map<int, Automaton> reference_populations = {
    {0, Automaton(0, 0.0, 0.0, 0.0, 0.0, 0.0)}, 
    {0x100, Automaton(0x100, b_growth, 0.0, 0.0, 0.0, 0.0)}, 
    {0x110, Automaton(0x110, b_growth, conjugation_010, s_loss, 0.0, 0.0)}, 
    {0x190, Automaton(0x190, b_growth, conjugation, s_loss, 0.0,p_resistance)}, 
    {0x900, Automaton(0x900, b_growth, 0.0, 0.0, b_resistance, 0.0)}, 
    {0x910, Automaton(0x910, b_growth, conjugation_010, s_loss, b_resistance, 0.0)}, 
    {0x990, Automaton(0x990, b_growth, conjugation, s_loss, b_resistance, p_resistance)},
    //lysogens without phage resistance 
    {0x101, Automaton(0x101, b_growth, 0.0, 0.0, 0.0, 0.0, burst_size_001, infection_rate, lysis_rate_001, lysogenization_rate_001, v_degradation_rate, v_loss_rate, 0.0)}, 
    {0x111, Automaton(0x111, b_growth, conjugation_010, s_loss, 0.0, 0.0, burst_size_001, infection_rate, lysis_rate_001, lysogenization_rate_001, v_degradation_rate, v_loss_rate, 0.0)}, 
    {0x191, Automaton(0x191, b_growth, conjugation, s_loss, 0.0,p_resistance, burst_size_001, infection_rate, lysis_rate_001, lysogenization_rate_001, v_degradation_rate, v_loss_rate, 0.0)}, 
    {0x901, Automaton(0x901, b_growth, 0.0, 0.0, b_resistance, 0.0, burst_size_001, infection_rate, lysis_rate_001, lysogenization_rate_001, v_degradation_rate, v_loss_rate, 0.0)}, 
    {0x911, Automaton(0x911, b_growth, conjugation_010, s_loss, b_resistance, 0.0, burst_size_001, infection_rate, lysis_rate_001, lysogenization_rate_001, v_degradation_rate, v_loss_rate, 0.0)}, 
    {0x991, Automaton(0x991, b_growth, conjugation, s_loss, b_resistance, p_resistance, burst_size_001, infection_rate, lysis_rate_001, lysogenization_rate_001, v_degradation_rate, v_loss_rate, 0.0)},
    //lysogens with phage resistance
    {0x109, Automaton(0x109, b_growth, 0.0, 0.0, 0.0, 0.0, burst_size, infection_rate, lysis_rate, lysogenization_rate, v_degradation_rate, v_loss_rate, v_resistance)}, 
    {0x119, Automaton(0x119, b_growth, conjugation_010, s_loss, 0.0, 0.0, burst_size, infection_rate, lysis_rate, lysogenization_rate, v_degradation_rate, v_loss_rate, v_resistance)}, 
    {0x199, Automaton(0x199, b_growth, conjugation, s_loss, 0.0,p_resistance, burst_size, infection_rate, lysis_rate, lysogenization_rate, v_degradation_rate, v_loss_rate, v_resistance)}, 
    {0x909, Automaton(0x909, b_growth, 0.0, 0.0, b_resistance, 0.0, burst_size, infection_rate, lysis_rate, lysogenization_rate, v_degradation_rate, v_loss_rate, v_resistance)}, 
    {0x919, Automaton(0x919, b_growth, conjugation_010, s_loss, b_resistance, 0.0, burst_size, infection_rate, lysis_rate, lysogenization_rate, v_degradation_rate, v_loss_rate, v_resistance)}, 
    {0x999, Automaton(0x999, b_growth, conjugation, s_loss, b_resistance, p_resistance, burst_size, infection_rate, lysis_rate, lysogenization_rate, v_degradation_rate, v_loss_rate, v_resistance)}
    }; 

    //Initialize the CA. Note that [0][col], [101][col], [row][0], [row][101] are the boundaries, whose states are usually fixed.
    for(unsigned row=1, p=0; row<rows_n+1; ++row){
      for(unsigned col=1; col<cols_n+1; ++col, ++p){
        ca_grid->cell(row,col) = reference_populations[populations[p]];
        if (!virions[p].empty()) virion_grid->cell(row,col) = virions[p];
      }
    }


    //BLOCK 3 - OUTPUT BLOCK
    std::string filename = experiment_tag + ".txt";
    std::ofstream final_counts(filename, std::ios::out);
    if ((!final_counts.is_open()) && output_data_flag) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    }
    unsigned relevant_pops[] = {0, 
    0x100, 0x110, 0x190, 0x900, 0x910, 0x990, 
    0x101, 0x111, 0x191, 0x901, 0x911, 0x991,
    0x109, 0x119, 0x199, 0x909, 0x919, 0x999,
    0x001, 0x009};
    final_counts << "time\t000CA\t100CA\t110CA\t190CA\t900CA\t910CA\t990CA\t101CA\t111CA\t191CA\t901CA\t911CA\t991CA\t109CA\t119CA\t199CA\t909CA\t919CA\t999CA\tvirions\tconjS\tconjR\tconjTot\tlysS\tlysR\tlysTot\tlysgS\tlysgR\tlysgTot\n";

  // set up state-to-color indices
  std::unordered_map<unsigned, unsigned> pixel_colors = {
    {0x000, 0},
    {0x100, 1}, {0x101, 7}, {0x109, 13}, 
    {0x110, 2}, {0x111, 8}, {0x119, 14},
    {0x190, 3}, {0x191, 9}, {0x199, 15}, 
    {0x900, 4}, {0x901, 10}, {0x909, 16}, 
    {0x910, 5}, {0x911, 11}, {0x919, 17}, 
    {0x990, 6}, {0x991, 12}, {0x999, 18}, 
    {0x001, 19}, {0x009, 20}, {0,0}, 
    //following types should not occur in simulation - used to detect errors in code (if these miscolored pixels are present, there is a problem)
    {0x105, 21},
    {0x115, 22}, 
    {0x195, 23}, 
    {0x905, 24}, 
    {0x915, 25}, 
    {0x995, 26},
    {0x10D, 27},  
    {0x11D, 28}, 
    {0x19D, 29},
    {0x90D, 30}, 
    {0x91D, 31},
    {0x99D, 32}
  };

  //OUTPUT CONFIG EARLY IN CASE RUN DOES NOT COMPLETE

  std::string config_name = experiment_tag + "_config.txt";
  std::ofstream config_file(config_name, std::ios::out);
  if ((!config_file.is_open()) && output_data_flag) {
  std::cerr << "Failed to open file: " << config_name << std::endl;
  }
  else
  {
    config_file << "SIMULATION OPTIONS\n";
    for (auto each:s_ops) config_file << each.first << ":\t" << each.second << "\n";
    config_file << "\nAUTOMATON PARAMS\n";
    for (auto each:a_prms) config_file << each.first << ":\t" << each.second << "\n";
    config_file << "\nBACTERIAL INFLUX RATES\n";
    for (auto each:bacteria_types) config_file << std::hex << each << ":\t" << influx_populations[pixel_colors[each]] << "\n";
    config_file.close();
  }

  //population tracking for outputs
    std::unordered_map<unsigned,unsigned> bacterial_counts;
    std::unordered_map<unsigned,unsigned> initial_bacterial_counts;

  /* Set parameters needed to display the CA. For this demo, we create
     only one panel */
  std::vector<CashPanelInfo> panel_info(2);

  //  Set a display panel size, set origin coordinate of display panel, and instantiate display object
  panel_info[0].n_row = rows_n;
  panel_info[0].n_col = cols_n;

  panel_info[0].o_row = 0;
  panel_info[0].o_col = 0;

  panel_info[1].n_row = rows_n;
  panel_info[1].n_col = cols_n;

  panel_info[1].o_row = 0;
  panel_info[1].o_col = cols_n;


  CashDisplay* display_p = nullptr;
  try{
    /* Window size is set to grid_dimension*grid_dimension. Only one window is allowed to open. */
    display_p = new CashDisplay(rows_n,cols_n*2,panel_info,scale);
  }catch(std::bad_alloc){
    std::cerr << "main(): Error, memory exhaustion" << std::endl;
    exit(-1);
  }

  /* If an X window needed, initialize things */
  if(show_display_flag){
    display_p->open_window("Demo CA");
  }
  
  /* If PNG slides are needed, initialize things */
  if(make_movie_flag){
    display_p->open_png(experiment_tag + "_movie");
  }

  /* Update the display */

    unsigned char color;
    unsigned char virion_color;
    // for(unsigned row=1; row<rows_n+1; ++row) {
    //   for(unsigned col=1; col<cols_n+1; ++col) {
    //     color = pixel_colors[ca_grid->cell(row,col).state];
    //     if (color == 29) testing(ca_grid->cell(row,col).state);
	  //     display_p->put_pixel(0,row,col,22);
    //   }
    // }

  //BLOCK 4 - Run the simulation

  //Instantiate a random number generator
  std::seed_seq seed({the_seed});
  std::mt19937_64 random(seed);
  std::uniform_int_distribution<int> uniform6(0,6);
  std::uniform_real_distribution<double> uniform1(0,1.0);
  std::uniform_int_distribution<int> die8(1,8);
  std::uniform_int_distribution<int> rand_grid(1,rows_n);
  std::uniform_int_distribution<int> rand_grid_sq(1,rows_n*cols_n);
  std::uniform_int_distribution<int> rand_virion(0,std::pow(2,31)); //choose a rand_virion from each cell

  //progress the model

  bool end_early_trigger = 0;

  //These averages are used to update the a_prms values at the end of the simulation
  double mean_evolved_lysis_tot = lysogenization_rate;
  double mean_evolved_lysogenization_tot = lysis_rate;
  double mean_evolved_conjugation_tot = conjugation;

  double mean_evolved_lysis_sus = lysis_rate_001;
  double mean_evolved_lysogenization_sus = lysogenization_rate_001;
  double mean_evolved_conjugation_sus = conjugation_010;

  double mean_evolved_lysis_res = lysis_rate;
  double mean_evolved_lysogenization_res = lysogenization_rate;
  double mean_evolved_conjugation_res = conjugation;


  for (unsigned time = 0; time<max_time; ++time) {
    if (time % 10000 == 0) std::cout << "TIME: " << time << std::endl;
    const bool draw_now_time = (time % draw_gap_time == 0) && show_display_flag;
    bool output_data_trigger = false;
    bool make_movie_trigger = false;
    bool init_movie_trigger = false; 
    bool micro_movie_trigger = false;

    if(make_movie_flag) {
      if(time < micro_movie_lim || time == micro_movie_lim) {
        micro_movie_trigger = (time % micro_movie_gap == 0);
      }
      if(time < init_movie_lim || time == init_movie_lim) {
        init_movie_trigger = (time % init_movie_gap == 0);
      }
      make_movie_trigger = (time % movie_gap == 0);
    }

    if(output_data_flag) {
      output_data_trigger = (time % output_gap == 0);
    }

    if (output_data_trigger) {
      for (auto pop: relevant_pops) {
        bacterial_counts[pop] = 0;
      }

      double cumulative_lysis_rate_sus = 0.0;
      double cumulative_lysis_rate_res = 0.0;
      double cumulative_lysogenization_rate_sus = 0.0;
      double cumulative_lysogenization_rate_res = 0.0;
      double cumulative_conjugation_rate_sus = 0.0;
      double cumulative_conjugation_rate_res = 0.0;
      double total_lysogens_sus = 0.0;
      double total_lysogens_res = 0.0;
      double total_plasmid_carriers_sus = 0.0;
      double total_plasmid_carriers_res = 0.0;
      for(unsigned row=1; row<rows_n+1; ++row) {
        for(unsigned col=1; col<cols_n+1; ++col) {
          bacterial_counts[ca_grid->cell(row,col).get_state()]++;
          bacterial_counts[0x001] += virion_grid->cell(row,col).size();
          if (ca_grid->cell(row,col).is_lysogen()) 
          {
            if (ca_grid->cell(row,col).has_resistant_phage()) 
            {
              cumulative_lysis_rate_res += ca_grid->cell(row,col).get_lysis_rate(false);
              cumulative_lysogenization_rate_res += ca_grid->cell(row,col).get_lysogenization_rate(false);
              total_lysogens_res += 1.0;
            }
            else
            {
              cumulative_lysis_rate_sus += ca_grid->cell(row,col).get_lysis_rate(false);
              cumulative_lysogenization_rate_sus += ca_grid->cell(row,col).get_lysogenization_rate(false);
              total_lysogens_sus += 1.0;
            }
          }
          if (ca_grid->cell(row,col).has_plasmid())
          {
            if (ca_grid->cell(row,col).has_resistant_plasmid()) 
            {
              cumulative_conjugation_rate_res += ca_grid->cell(row,col).get_plasmid_conjugation_rate(false);
              total_plasmid_carriers_res += 1.0;
            }
            else
            {
              cumulative_conjugation_rate_sus += ca_grid->cell(row,col).get_plasmid_conjugation_rate(false);
              total_plasmid_carriers_sus += 1.0;
            }
          }
        }
      }

      mean_evolved_lysis_sus = cumulative_lysis_rate_sus / std::max(1.0,total_lysogens_sus); 
      mean_evolved_lysis_res = cumulative_lysis_rate_res / std::max(1.0,total_lysogens_res);
      mean_evolved_lysis_tot = (cumulative_lysis_rate_sus + cumulative_lysis_rate_res) / std::max(1.0, total_lysogens_sus + total_lysogens_res);
      
      mean_evolved_lysogenization_sus = cumulative_lysogenization_rate_sus / std::max(1.0,total_lysogens_sus);
      mean_evolved_lysogenization_res = cumulative_lysogenization_rate_res / std::max(1.0,total_lysogens_res);
      mean_evolved_lysogenization_tot = (cumulative_lysogenization_rate_sus + cumulative_lysogenization_rate_res) / std::max(1.0, total_lysogens_sus + total_lysogens_res);
      
      mean_evolved_conjugation_sus = cumulative_conjugation_rate_sus / std::max(1.0,total_plasmid_carriers_sus);
      mean_evolved_conjugation_res = cumulative_conjugation_rate_res / std::max(1.0,total_plasmid_carriers_res);
      mean_evolved_conjugation_tot = (cumulative_conjugation_rate_sus + cumulative_conjugation_rate_res) / std::max(1.0, total_plasmid_carriers_sus + total_plasmid_carriers_res);

      final_counts << time;
      for (auto pop: relevant_pops) {
        if (pop != 0x009) final_counts << "\t" << bacterial_counts[pop];
      }
      final_counts << "\t" << mean_evolved_conjugation_sus << "\t" << mean_evolved_conjugation_res << "\t" << mean_evolved_conjugation_tot;
      final_counts << "\t" << mean_evolved_lysis_sus << "\t" << mean_evolved_lysis_res << "\t" << mean_evolved_lysis_tot;
      final_counts << "\t" << mean_evolved_lysogenization_sus << "\t" << mean_evolved_lysogenization_res << "\t" << mean_evolved_lysogenization_tot;
      final_counts << std::endl;
    }

    /* If PNG slides are needed, draw things */
    if(make_movie_trigger || init_movie_trigger || micro_movie_trigger)
    {
      for(unsigned row=1; row<rows_n+1; ++row) {
        for(unsigned col=1; col<cols_n+1; ++col) {
          //output bacteria grid to movie (using color mapping for cell states)
          color = pixel_colors[ca_grid->cell(row,col).get_state()];
          display_p->put_pixel(0,row,col,color);
          //output virion grid to movie (using intensity scaling for virion count)
          virion_color = virion_grid->cell(row,col).size() + 40;
          if (virion_color > 43) virion_color = 43;  // Cap maximum intensity
          display_p->put_pixel(1,row,col,virion_color);
        }
      }
      display_p->draw_png();
    }
    //within each timestep, have gridsize opportunities to visit all pixels on the grid with equal probability of visiting any pixel per pixel visit
    for(unsigned pixel_visit =0; pixel_visit<gridsize; pixel_visit++) {

    bool show_display_trigger = (pixel_visit % draw_gap == 0) * show_display_flag * draw_now_time;
    antibiotic_trigger = (time % antibiotic_time_periodicity < antibiotic_time_limit) && (pixel_visit % antibiotic_pixel_periodicity < antibiotic_pixel_limit); 

    /* Update display */
    if (show_display_trigger) {
      for(unsigned row=1; row<rows_n+1; ++row) {
        for(unsigned col=1; col<cols_n+1; ++col) {
          color = pixel_colors[ca_grid->cell(row,col).get_state()];
          if ((ca_grid->cell(row,col).get_state() & 0x005) == 0x005) testing(ca_grid->cell(row,col).get_state());
          display_p->put_pixel(0,row,col,color);

          virion_color = virion_grid->cell(row,col).size() + 40;
          if (virion_color > 43) virion_color = 43;
          display_p->put_pixel(1,row,col,virion_color);
        }
      } 
    }

    /* If an X window needed, draw things */
    if(show_display_trigger)
    {
      display_p->draw_window();
        // delay program to see CA in real pixel_visit   
      if(test_flag) std::this_thread::sleep_for(std::chrono::milliseconds(delay_period));   
    }

     // this block controls whether to end simulation early. If simuation_opts["end_early"] == 1 AND either resistant chromosomes or susceptible chromosomes have gone extenct, then the simulation will end early
        if (end_early && (time > 0)) {
          if (((bacterial_counts[0x100] == 0) && (bacterial_counts[0x110] == 0) && (bacterial_counts[0x190] == 0)) || ((bacterial_counts[0x900] == 0) && (bacterial_counts[0x910] == 0) && (bacterial_counts[0x990] == 0))) {
            if (output_data_trigger && make_movie_trigger) end_early_trigger = 1;
          }
        }

    for (unsigned grid_i = 0; grid_i < 2; ++grid_i) {
      // choose a random cell and random neighbour
      unsigned row_i = rand_grid(random); 
      unsigned col_i = rand_grid(random);
      unsigned rand_nei = die8(random);

      double event_chance = uniform1(random);
      
        //initialize aliases for current cell and current neighbour - easier to read
      
      if (diffuse_switch) std::swap(ca_grid->cell(row_i,col_i),ca_grid->neigh_wrap(row_i,col_i,die8(random)));
      Automaton* cur_cell = &ca_grid->cell(row_i,col_i);
      Automaton* cur_nei = &ca_grid->neigh_wrap(row_i,col_i, rand_nei); //random mixing has been hard-coded - needs to be ammended later
      std::vector<Automaton>* cur_virions = &virion_grid->cell(row_i, col_i);
      
      if (random_mixing) 
      {
        cur_nei = &ca_grid->cell(rand_grid(random), rand_grid(random));
        cur_cell = &ca_grid->cell(rand_grid(random),rand_grid(random));
      }



      //define event probabilities and conditions for ease of reading the control blocks. Use lambda so that these conditions "automatically" update, in case they are used multiple times in a single loop iteration

      auto virion_present = [cur_virions]() -> bool {
        return !(cur_virions->empty());
      };

      auto bacterial_death_threshold = [&cur_cell, antibiotic_trigger,antibiotic_pressure,b_death]() -> double {
        return std::max(b_death + antibiotic_trigger*(antibiotic_pressure - cur_cell->get_plasmid_resistance() - cur_cell->get_bacterial_resistance() - cur_cell->get_phage_resistance()), b_death); 
      };

      auto conjugation_threshold = [&cur_cell, bacterial_death_threshold]() -> double {
        return cur_cell->get_plasmid_conjugation_rate() + bacterial_death_threshold();
      };

      auto lysis_threshold = [&cur_cell, conjugation_threshold]() -> double {
        return cur_cell->get_lysis_rate() + conjugation_threshold();
      };

      auto reproduction_threshold = [&cur_nei, r_cost, p_cost, b_growth]() -> double {
        return std::max(b_growth - (cur_nei->has_plasmid()  * p_cost) - (cur_nei->has_resistant_chromosome() * r_cost) - (cur_nei->has_resistant_plasmid() * r_cost) -
          (cur_nei->has_resistant_phage() * r_cost), 0.0);
      };

      auto plasmid_loss_threshold = [&cur_nei]() -> double {
        return cur_nei->get_plasmid_loss_rate();
      };

      auto phage_loss_threshold = [&cur_nei, plasmid_loss_threshold]() -> double {
        return plasmid_loss_threshold() + cur_nei->get_phage_loss_rate();
      };

      auto v_degradation_threshold = [&cur_virions, v_diffusion](int i) -> double {
        return (*cur_virions)[i].get_virion_degradation_rate() + v_diffusion;
      };

      auto v_infection_threshold = [&cur_virions, v_diffusion, v_degradation_threshold](int i) -> double {
        return (*cur_virions)[i].get_infection_rate() + v_degradation_threshold(i);
      };
      //TEMPORARY
      auto mutate_lysogenization = [delta_t] (Automaton* cell) -> void { 
        double temp_lysogenization = cell->get_lysogenization_rate(false);
        double mutated_lysogenization = cell->mutate(temp_lysogenization, 0.005*delta_t, 1.0);
        cell->set_lysogenization_rate(mutated_lysogenization);
      };

      //TEMPORARY
      auto mutate_background_lysis = [delta_t] (Automaton* cell) -> void { 
        double temp_lysis_rate = cell->get_lysis_rate(false);
        double mutated_lysogenization = cell->mutate(temp_lysis_rate, 0.005*delta_t, 0.2);
        cell->set_lysis_rate(mutated_lysogenization);
      };

      //TEMPORARY
      auto mutate_conjugation = [delta_t] (Automaton* cell) -> void { 
        double temp_conjugation_rate = cell->get_plasmid_conjugation_rate(false);
        double mutated_conjugation = cell->mutate(temp_conjugation_rate, 0.005*delta_t, 0.2);
        cell->set_plasmid_conjugation_rate(mutated_conjugation);
      };


      //v_infection_threshold, v_degradation_threshold, and v_diffusion must together be 1 or less


      //All virions present get an opportunity to infect current cell. If current cell is already infected, virion is simply wasted and disappears
      if (grid_i == 0) //handles virion gridpoints and bacterial gridpoints separately, to improve randomness of the the model
      {
        if (virion_present()) 
        {
          // random shuffle virions
          std::shuffle(cur_virions->begin(),cur_virions->end(), random);
          {
            //iterate through virions, giving each one a chance to infect cell, fail to infect, or degrade
            for (unsigned i = cur_virions->size(); i-- > 0; ) 
            {
              if (i >= cur_virions->size()) continue; // Ensure index is within bounds
              double virion_event_chance = uniform1(random);
              if (random_mixing) cur_cell = &ca_grid->cell(rand_grid(random),rand_grid(random));

              if (virion_event_chance < v_diffusion) 
              { 
                //diffuse virions by removing them from current cells virion vector and adding them to random neighbouring cell
                if (random_mixing) virion_grid->cell(rand_grid(random),rand_grid(random)).push_back(cur_virions->back());
                else virion_grid->neigh_wrap(row_i,col_i,die8(random)).push_back(cur_virions->back());

                cur_virions->pop_back();
              } 
              else if (virion_event_chance < v_degradation_threshold(i))
              {
                //virion degrades and ceases to exist
                cur_virions->pop_back();
              } 
              else if (virion_event_chance < v_infection_threshold(i))
              {
                //virion "tries to infect"
                if (cur_cell->is_alive())  
                {
              //actual attempt to infect only occurs if there is a living bacteria present
              if (!cur_cell->is_lysogen()) 
              {
                  //infection can only succeed if cell is not a lysogen
                  (*cur_virions)[i].viral_infect(cur_cell);
                  if (uniform1(random) >= cur_cell->get_lysogenization_rate())
                  {
                for (unsigned j = 0; j < cur_cell->get_burst_size(); j++) 
                {
                    cur_virions->push_back(cur_cell->make_virion());
                }
                cur_cell->die();
                  } 
                  cur_virions->pop_back();
              } else {
                  //if the cell is a lysogen, the virion is wasted while attempting to infect the cell, and ceases to exist
                  cur_virions->pop_back();
              }
                }
              } 
            }
          }
        }
      }
      else if (grid_i == 1)//bacterial gridpoints are chosen separately from virion gridpoints to improve randomness of the model
      {
          //Is bacteria alive? 
          if (cur_cell->is_alive()) {

            // bacteria can die spontaneously or due to antibiotics...
            if (event_chance < bacterial_death_threshold()) {
              if (death_lysis_switch) //lysis in response to death
              {
                for (unsigned i = 0; i < cur_cell->get_burst_size(); i++) 
                {
                  cur_virions->push_back(cur_cell->make_virion());
                  if (evolve_parameters_switch) 
                  {
                    mutate_background_lysis(&cur_virions->back());
                    mutate_lysogenization(&cur_virions->back());
                  }
                }
              }
              cur_cell->die();
            // ...plasmid can conjugate to a random neighbour...
            } else if (event_chance < conjugation_threshold()) 
            {
              if (cur_cell->has_plasmid() && cur_nei->is_alive() && (!cur_nei->has_plasmid())) 
              {
                cur_cell->conjugate(cur_nei);
                if (evolve_parameters_switch) mutate_conjugation(cur_nei); //TEMPORARY -- we considider, as a simplification, that conjugation rate mutates upon conjugation because plasmids have regulated copy numbers to some target number - if plasmid count is decreased within a cell beneath the target, a new plasmid should be created, which gives the new plasmid an opportunity to gain mutations during the plasmid creation process. In our model, although individual plasmids within a cell are not modelled explicitly, if a plasmid is passed to another cell, it means implicitly that the plasmid already exists, and has already undergone the plasmid copying process, and may have mutated already. At this point, we may consider whether indeed such a plasmid has mutations. For simplicity, we are assuming that mutated plasmids do not matter in the host in which they are first generated, since their copy number will at first be low; but they matter a lot in the new host which they infect, since they will be the founding template plasmid in the new host, from which new plasmids are copied.  
              }
            } else if (event_chance < lysis_threshold()) //lyse
            { 
              for (unsigned i = 0; i < cur_cell->get_burst_size(); i++) 
              {
                  cur_virions->push_back(cur_cell->make_virion());
                  if (evolve_parameters_switch) 
                  {
                    mutate_background_lysis(&cur_virions->back());
                    mutate_lysogenization(&cur_virions->back());
                  }
              }
              cur_cell->die();
            }
          } else if (cur_cell->is_empty()) {
            
            //...empty cells allow living neighbours to reproduce onto them
            if (cur_nei->is_alive() && (event_chance < reproduction_threshold())) {
              double loss_chance = uniform1(random);
              if ((loss_chance > phage_loss_threshold()) || (loss_chance == phage_loss_threshold())) 
              {
                cur_cell->get_offspring(cur_nei); //reproduce bacteria
              }
              else if (loss_chance < plasmid_loss_threshold())
              {
                cur_cell->get_offspring(cur_nei,true); //reproduce bacteria, lose plasmid
              }        
              else 
              {
                cur_cell->get_offspring(cur_nei, false, true); //lose phage
              }
              //if current cell is plasmid-bearing, evolve conjugation - TEMPORARY
              if (cur_cell->has_plasmid() && evolve_parameters_switch)
              {
                mutate_conjugation(cur_cell);
              }
              //if current cell is lysogen, evolve ab_lysis_response - TEMPORARY
              if (cur_cell->is_lysogen() && evolve_parameters_switch)
              {
                mutate_background_lysis(cur_cell);
                mutate_lysogenization(cur_cell);
              }

              // ...or empty cells can be seeded onto by constant bacterial influx
            } else if (event_chance < (influx_rate + reproduction_threshold()) && event_chance > reproduction_threshold()) 
            {
              influx_chance = reproduction_threshold();
              for (unsigned pop = 0; pop < influx_populations.size(); pop++) {
                influx_chance += influx_populations[pop];
                if (event_chance <  influx_chance) {
              cur_cell->get_offspring(&reference_populations[bacteria_types[pop]]);
              break;
                }
              }
              if (evolve_parameters_switch)
              {
                // To newly influxed cells, apply rates similar to average evolved rates - this gives invading bacteria comparable fitness and a better chance of invading. In this way, if invasion fails, it is not just because the invading bacteria have less fit evolved parameters than the established population
                if (cur_cell->is_lysogen()) 
               if (cur_cell->is_lysogen()) 
                {
                  if (mean_evolved_lysogenization_tot > 0.0)  cur_cell->set_lysogenization_rate(mean_evolved_lysogenization_tot);
                  else cur_cell->set_lysogenization_rate(lysogenization_rate_t0);
                  if (mean_evolved_lysis_tot > 0.0) cur_cell->set_lysis_rate(mean_evolved_lysis_tot);
                  else cur_cell->set_lysis_rate(lysis_rate_t0);

                  mutate_background_lysis(cur_cell);
                  mutate_lysogenization(cur_cell);

                } 
                if (cur_cell->has_plasmid()) {
                  if (mean_evolved_conjugation_tot > 0.0) cur_cell->set_plasmid_conjugation_rate(mean_evolved_conjugation_tot);
                  else cur_cell->set_plasmid_conjugation_rate(conjugation_t0);
                  mutate_conjugation(cur_cell);
                }
              }
            }
          }
        }
      }
    }
      if (end_early_trigger) break;
      }
      if (make_movie_flag || output_data_flag) {
        final_counts << "SIMULATION OPTIONS\n";
        for (auto each:s_ops) final_counts << each.first << ":\t" << each.second << "\n";
        final_counts << "\nAUTOMATON PARAMS\n";
        for (auto each:a_prms) final_counts << each.first << ":\t" << each.second << "\n";
        final_counts << "\nBACTERIAL INFLUX RATES\n";
        for (auto each:bacteria_types) final_counts << std::hex << each << ":\t" << influx_populations[pixel_colors[each]] << "\n";
  }
  final_counts.close();

  //take a snapshot of final state to return. In this way, the result of testruns can be used as the starting point for subsequent runs
  std::vector<unsigned> final_state(gridsize);
  for(unsigned row=1, i=0; row<rows_n+1; ++row) {
    for(unsigned col=1; col<cols_n+1; ++col, ++i) {
      final_state[i] = ca_grid->cell(row,col).get_state();
    }
  }

  delete ca_grid;
  delete display_p;
  delete virion_grid;
  std::cout << "DONE!\n";

  //updating evolved rates in a_prms to pass to the next simulation
  a_prms["lysis_rate"] = mean_evolved_lysis_res;
  a_prms["lysogenization_rate"] = mean_evolved_lysogenization_res;
  a_prms["conjugation"] = mean_evolved_conjugation_res;

  a_prms["lysis_rate_001"] = mean_evolved_lysis_sus;
  a_prms["lysogenization_rate_001"] = mean_evolved_lysogenization_sus;
  a_prms["conjugation_010"] = mean_evolved_conjugation_sus;

std::cout << "cj " << mean_evolved_conjugation_sus << " lys " << " " << mean_evolved_lysis_sus << " lysg " << mean_evolved_lysogenization_sus << std::endl;
  
  return final_state;
}
