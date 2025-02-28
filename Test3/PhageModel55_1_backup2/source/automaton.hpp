#include <iostream>
#include <vector>
#include <array>

#ifndef AUTOMATON
#define AUTOMATON

class Automaton {
private:
  //state: 0 = nothing. hex hundreds is bacterial column; hex tens is plasmid column; hex ones is phage column; for each hex column, there are four binary digits, of which last (right-most) binary digit is "existence" or "life" flag, and first binary digit is flag which shows whether that genetic element carries antibiotic resistence. EXAMPLE: 0x190 means a living bacteria with a non-resistant chromosome (described by the first digit, "1", or in binary "0001") which carries a resistant plasmid (described by "9", or "1001" in binary) and no phage (described by 0)
  static double delta_t;
  unsigned state; 

  double bacterial_growth_rate; 
  double bacterial_resistance;

  double plasmid_conjugation_rate;
  double plasmid_loss_rate;
  double plasmid_resistance;

  //Phage parameters
  double phage_resistance;
  unsigned burst_size;
  
  double infection_rate;
  double lysis;
  double lysogenization;
  double virion_degradation;
  double phage_loss;
  double ab_lysis_response; 
public:
  
  // double phage_growth_rate;
  // double phage_lysis_rate;
  // int phage_virions_n;
  // double virion_persistance; 

  Automaton();
  //bacterial initialization
  Automaton(unsigned, double, double, double, double, double); 
  //initialize lysogen
  Automaton (unsigned state_in, double growth, double conjugation, double loss, double b_resistance, double p_resistance, 
  unsigned burst_size_in, 
  double infection_rate_in, 
  double lysis_in, 
  double lysogenization_in, 
  double virion_degradation_in,
  double phage_loss_in, 
  double phage_resistance_in);

    Automaton (unsigned state_in, double growth, double conjugation, double loss, double b_resistance, double p_resistance, 
  unsigned burst_size_in, 
  double infection_rate_in, 
  double lysis_in, 
  double lysogenization_in, 
  double virion_degradation_in,
  double phage_loss_in, 
  double phage_resistance_in,
  double ab_lysis_response_in);
  
  //virion initialization 
  int init_virion(unsigned state_in, unsigned burst_size_in, double lysogenization_in, double lysis_in, double virion_degradation_in, double infection_rate_in, double phage_loss_in); 
  int init_virion(unsigned state_in, unsigned burst_size_in, double lysogenization_in, double lysis_in, double virion_degradation_in, double infection_rate_in, double phage_loss_in, double phage_resistance_in); 
  int init_virion(unsigned state_in, unsigned burst_size_in, double lysogenization_in, double lysis_in, double virion_degradation_in, double infection_rate_in, double phage_loss_in, double phage_resistance_in, double ab_lysis_response_in); 

  Automaton make_virion();
  double static mutate(double mutate_value, double mutate_limit, double max_mutate_value);

  int set(unsigned, double, double, double, double, double);
  int die();
  int get_phage(); //no origin of phage is specified, so phage is lost
  int get_phage(Automaton* neighbour);
  int reproduce_phage(Automaton* neighbour);
  int viral_infect(Automaton* neighbour);

  int conjugate(Automaton* neighbour);
  int reproduce(Automaton* neighbour);
  int reproduce_lose(Automaton* neighbour);
  int get_offspring(Automaton* neighbour);
  int get_offspring(Automaton* neighbour, bool lose_plasmid, bool lose_phage = false);
  int get_offspring_lose(Automaton* neighbour);

  bool is_alive();
  bool is_empty();
  bool has_plasmid();
  bool is_lysogen();
  bool has_resistant_plasmid();
  bool has_resistant_chromosome();
  bool has_resistant_phage();

  //get data about bacteria
  unsigned get_state(); 
  double get_bacterial_growth_rate(); 
  double get_bacterial_resistance();

  //get data about plasmid
  double get_plasmid_conjugation_rate(bool use_delta_t = true);
  double get_plasmid_loss_rate();
  double get_phage_loss_rate();
  double get_plasmid_resistance();

  //get data about phage
  double get_phage_resistance();
  unsigned get_burst_size();
  
  double get_infection_rate();
  double get_lysis_rate(bool use_delta_t = true);
  double get_lysogenization_rate(bool use_delta_t = true);
  double get_virion_degradation_rate();
  double get_ab_lysis_response(bool delta_t = true);

    //set data about bacteria
  void set_state(unsigned); 
  void set_bacterial_growth_rate(double); 
  void set_bacterial_resistance(double);

  //set data about plasmid
  void set_plasmid_conjugation_rate(double);
  void set_plasmid_loss_rate(double);
  void set_plasmid_resistance(double);

  //set data about phage
  void set_phage_resistance(double);
  void set_burst_size(unsigned);
  
  void set_infection_rate(double);
  void set_lysis_rate(double);
  void set_lysogenization_rate(double);
  void set_virion_degradation_rate(double);
  void set_ab_lysis_response(double);

  static double get_delta_t()
  {
    return delta_t;
  }

  static void set_delta_t(double value)
  {
    delta_t = value;
  }

};

#endif

