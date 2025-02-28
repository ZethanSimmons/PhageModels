#include <vector>
#include "automaton.hpp"
#include <random>
#include <cmath>

double Automaton::delta_t = 1.0;

Automaton::Automaton () {
    state = 0;
    bacterial_growth_rate = 0.0;
    plasmid_conjugation_rate = 0.0;
    plasmid_loss_rate = 0.0;
    bacterial_resistance = 0.0;
    plasmid_resistance = 0.0;

    //phage parameters
    phage_resistance = 0.0;

    burst_size = 0;
    lysis = 0.0;
    lysogenization = 0.0;
    infection_rate = 0.0;
    virion_degradation = 0.0;
    phage_loss = 0.0;
    ab_lysis_response = 0.0;
}

//initialize a bacteria
Automaton::Automaton (unsigned state_in, double growth, double conjugation, double loss, double b_resistance, double p_resistance) {
    state = state_in;
    bacterial_growth_rate = growth;
    plasmid_conjugation_rate = conjugation;
    plasmid_loss_rate = loss;
    bacterial_resistance = b_resistance;
    plasmid_resistance = p_resistance;

    //phage parameters
    phage_resistance = 0.0;

    burst_size = 0;
    infection_rate = 0.0;
    lysis = 0.0;
    lysogenization = 0.0;
    virion_degradation = 0.0;
    phage_loss = 0.0;
    ab_lysis_response = 0.0;
}

Automaton::Automaton (unsigned state_in, double growth, double conjugation, double loss, double b_resistance, double p_resistance, 
unsigned burst_size_in, 
double infection_rate_in, 
double lysis_in, 
double lysogenization_in, 
double virion_degradation_in, 
double phage_loss_in,
double phage_resistance_in) {
    state = state_in;
    bacterial_growth_rate = growth;
    plasmid_conjugation_rate = conjugation;
    plasmid_loss_rate = loss;
    bacterial_resistance = b_resistance;
    plasmid_resistance = p_resistance;

    //phage parameters
    phage_resistance = phage_resistance_in;

    burst_size = burst_size_in;
    infection_rate = infection_rate_in;
    lysis = lysis_in;
    lysogenization = lysogenization_in;
    virion_degradation = virion_degradation_in;
    phage_loss = phage_loss_in;
}

Automaton::Automaton (unsigned state_in, double growth, double conjugation, double loss, double b_resistance, double p_resistance, 
unsigned burst_size_in, 
double infection_rate_in, 
double lysis_in, 
double lysogenization_in, 
double virion_degradation_in, 
double phage_loss_in,
double phage_resistance_in,
double ab_lysis_response_in) 
{
    Automaton(state_in, growth, conjugation, loss, b_resistance,p_resistance, 
        burst_size_in, 
        infection_rate_in, 
        lysis_in, 
        lysogenization_in, 
        virion_degradation_in, 
        phage_loss_in,
        phage_resistance_in);
    ab_lysis_response = ab_lysis_response_in;
}

//initialize virion

int Automaton::init_virion (unsigned state_in, unsigned burst_size_in, double lysogenization_in, double lysis_in, double virion_degradation_in, double infection_rate_in, double phage_loss_in) {
    state = state_in;
    burst_size = burst_size_in;
    lysogenization = lysogenization_in;
    lysis = lysis_in;
    virion_degradation = virion_degradation_in;
    infection_rate = infection_rate_in;
    phage_loss = phage_loss_in;

    phage_resistance = 0.0;

    //bacterial parameters
    bacterial_growth_rate = 0.0;
    plasmid_conjugation_rate = 0.0;
    plasmid_loss_rate = 0.0;
    bacterial_resistance = 0.0;
    plasmid_resistance = 0.0;

    return 0;
    
}

int Automaton::init_virion (unsigned state_in, unsigned burst_size_in, double lysogenization_in, double lysis_in, double virion_degradation_in, double infection_rate_in, double phage_loss_in, double phage_resistance_in) {
    state = state_in;
    burst_size = burst_size_in;
    lysogenization = lysogenization_in;
    lysis = lysis_in;
    virion_degradation = virion_degradation_in;
    infection_rate = infection_rate_in;
    phage_loss = phage_loss_in;

    phage_resistance = phage_resistance_in;

    //bacterial parameters
    bacterial_growth_rate = 0.0;
    plasmid_conjugation_rate = 0.0;
    plasmid_loss_rate = 0.0;
    bacterial_resistance = 0.0;
    plasmid_resistance = 0.0;

    return 0;
    
}

int Automaton::init_virion (unsigned state_in, unsigned burst_size_in, double lysogenization_in, double lysis_in, double virion_degradation_in, double infection_rate_in, double phage_loss_in, double phage_resistance_in, double ab_lysis_response_in) {

    init_virion (state_in, burst_size_in, lysogenization_in,  lysis_in, virion_degradation_in, infection_rate_in, phage_loss_in, phage_resistance_in);

    ab_lysis_response = ab_lysis_response_in;

    return 0;
    
}


int Automaton::get_phage() {
    //with no origin of phage specified, phage is lost
    state = (0xFF0 & state);
    burst_size = 0;
    lysogenization = 0;
    lysis = 0;
    virion_degradation = 0;
    infection_rate = 0;
    phage_resistance = 0;
    phage_loss = 0;

    return 0;
}

int Automaton::get_phage(Automaton* neighbour) {
    state += (0x00F & neighbour->state);
    burst_size = neighbour->burst_size;
    lysogenization = neighbour->lysogenization;
    lysis = neighbour->lysis;
    virion_degradation = neighbour->virion_degradation;
    infection_rate = neighbour->infection_rate;
    phage_resistance = neighbour->phage_resistance;
    phage_loss = neighbour->phage_loss;
    ab_lysis_response = neighbour->ab_lysis_response;

    return 0;
}

int Automaton::reproduce_phage(Automaton* neighbour) {
    neighbour->state += (0x00F & state);
    neighbour->burst_size = burst_size;
    neighbour->lysogenization = lysogenization;
    neighbour->lysis = lysis;
    neighbour->virion_degradation = virion_degradation;
    neighbour->infection_rate = infection_rate;
    neighbour->phage_resistance = phage_resistance;
    neighbour->phage_loss = phage_loss;
    neighbour->ab_lysis_response = ab_lysis_response;

    return 0;
}

int Automaton::viral_infect(Automaton* neighbour) {
    //infect cell if not already infected
    reproduce_phage(neighbour);
    die();
    return 0;
}

Automaton Automaton::make_virion() {
    Automaton new_virion{};
    new_virion.init_virion((state & 0x00F), burst_size, lysogenization, lysis, virion_degradation, infection_rate, phage_loss, phage_resistance);
    return new_virion;
}



double Automaton::mutate(double mutate_value, double mutate_limit, double max_mutate_value) {
    // Create a random number generator
    static std::random_device rd;
    static std::mt19937 gen(rd());
    
    // Create a uniform distribution for the mutation value
    std::uniform_real_distribution<> dis(0, mutate_limit);
    
    // Generate a random mutation
    double mutation = dis(gen);
    
    // 50% chance of being negative
    if (std::bernoulli_distribution(0.5)(gen)) {
        mutation = -mutation;
    }
    
    // Calculate the new value
    double new_value = mutate_value + mutation;
    
    // Handle reversal for values outside the allowed range
    if (new_value < 0) {
        new_value = std::abs(new_value);  // Reverse from negative
    } else if (new_value > max_mutate_value) {
        double excess = new_value - max_mutate_value;
        new_value = max_mutate_value - excess;  // Reverse from top
    }
    
    // In case of double reversal, ensure we're within bounds
    return std::max(0.0, std::min(new_value, max_mutate_value));
}

int Automaton::set(unsigned state_in = 0, double growth = 0, double conjugation = 0, double loss = 0, double b_resistance = 0, double p_resistance = 0) {
    state = state_in;
    bacterial_growth_rate = growth;
    plasmid_conjugation_rate = conjugation;
    plasmid_loss_rate = loss;
    bacterial_resistance = b_resistance;
    plasmid_resistance = p_resistance;
    return 0;
}

int Automaton::die()
{
    state = 0x0;
    bacterial_growth_rate = 0;
    bacterial_resistance = 0;
    plasmid_conjugation_rate = 0;
    plasmid_loss_rate = 0;
    plasmid_resistance = 0;

    //phage parameters
    burst_size = 0;
    lysis = 0.0;
    lysogenization = 0.0;
    virion_degradation = 0.0;
    
    return 0;
}

int Automaton::conjugate(Automaton* neighbour) 
{
    neighbour->state = (neighbour->state) | (0x0F0 & state);
    neighbour->plasmid_conjugation_rate = plasmid_conjugation_rate;
    neighbour->plasmid_loss_rate = plasmid_loss_rate;
    neighbour->plasmid_resistance = plasmid_resistance;  
    return 0;
}

 int Automaton::reproduce(Automaton* neighbour)
  {
    neighbour->state = (0xFF0 & state);
    neighbour->bacterial_growth_rate = bacterial_growth_rate;
    neighbour->bacterial_resistance = bacterial_resistance;
    neighbour->plasmid_loss_rate = plasmid_loss_rate;
    neighbour->plasmid_conjugation_rate = plasmid_conjugation_rate;
    neighbour->plasmid_resistance = plasmid_resistance;

    //phage parameters
    reproduce_phage(neighbour);

    return 0;
  }

  int Automaton::reproduce_lose(Automaton* neighbour)
  {
    //give bacterial aspects but nothing else
    neighbour->state = (0xF00 & state);
    neighbour->bacterial_growth_rate = bacterial_growth_rate;
    neighbour->bacterial_resistance = bacterial_resistance;
    neighbour->plasmid_loss_rate = 0;
    neighbour->plasmid_conjugation_rate = 0;
    neighbour->plasmid_resistance = 0;

    //phage parameters
    reproduce_phage(neighbour);

    return 0;
  }

int Automaton::get_offspring(Automaton* neighbour)
{
    state = (0xFF0 & neighbour->state); //phage state is handled by get_phage() function
    bacterial_growth_rate = neighbour->bacterial_growth_rate;
    bacterial_resistance = neighbour->bacterial_resistance;
    plasmid_conjugation_rate = neighbour->plasmid_conjugation_rate;
    plasmid_loss_rate = neighbour->plasmid_loss_rate;
    plasmid_resistance = neighbour->plasmid_resistance;

    //phage parameters
    get_phage(neighbour);

    return 0;
}

int Automaton::get_offspring(Automaton* neighbour, bool lose_plasmid, bool lose_phage)
{
    state = (0xFF0 & neighbour->get_state()); 
    bacterial_growth_rate = neighbour->bacterial_growth_rate;
    bacterial_resistance = neighbour->bacterial_resistance;

    if (lose_plasmid)
    {
        state = (state & 0xF00); //lose plasmid
        plasmid_conjugation_rate = 0;
        plasmid_loss_rate = 0;
        plasmid_resistance = 0;
    } 
    else
    {
        plasmid_conjugation_rate = neighbour->plasmid_conjugation_rate;
        plasmid_loss_rate = neighbour->plasmid_loss_rate;
        plasmid_resistance = neighbour->plasmid_resistance;
    }
    //phage parameters
    if (lose_phage) get_phage();
    else get_phage(neighbour);
    

    return 0;
}

int Automaton::get_offspring_lose(Automaton* neighbour)
{
    //get bacterial aspects of neighbour but lose plasmid
    state = (0xF00 & neighbour->state);
    bacterial_growth_rate = neighbour->bacterial_growth_rate;
    bacterial_resistance = neighbour->bacterial_resistance;
    plasmid_conjugation_rate = 0;
    plasmid_loss_rate = 0;
    plasmid_resistance = 0;

    //phage parameters
    get_phage(neighbour);

    return 0;
}

bool Automaton::is_alive()
{
    return (state & 0x100) == 0x100;
}

bool Automaton::is_empty()
{
    return (state & 0x100) == 0x000;
}

bool Automaton::has_plasmid()
{
    return (state & 0x010) == 0x010;
}

bool Automaton::is_lysogen()
{
    return (state & 0x001) == 0x001;
}

bool Automaton::has_resistant_plasmid()
{
    return (state & 0x090) == 0x090;
}

bool Automaton::has_resistant_chromosome()
{
    return (state & 0x900) == 0x900;
}

bool Automaton::has_resistant_phage()
{
    return (state & 0x009) == 0x009;
}

 //get data about bacteria
unsigned Automaton::get_state()
{
    return state;
} 

double Automaton::get_bacterial_growth_rate()
{
    return get_delta_t() * bacterial_growth_rate;
}

double Automaton::get_bacterial_resistance()
{
    return get_delta_t() * bacterial_resistance;
}

//get data about plasmid
double Automaton::get_plasmid_conjugation_rate(bool use_delta_t)
{
    if (use_delta_t)
    {
    return get_delta_t() * plasmid_conjugation_rate;}
    else
    {
        return plasmid_conjugation_rate;
    }
}

double Automaton::get_plasmid_loss_rate()
{
    return get_delta_t() * plasmid_loss_rate;
}

double Automaton::get_phage_loss_rate()
{
    return get_delta_t() * phage_loss;
}

double Automaton::get_plasmid_resistance()
{
    return get_delta_t() * plasmid_resistance;
}

//get data about phage
double Automaton::get_phage_resistance()
{
    return get_delta_t() * phage_resistance;
}

unsigned Automaton::get_burst_size()
{
    return burst_size;
}

double Automaton::get_infection_rate()
{
    return get_delta_t() * infection_rate;
}

double Automaton::get_lysis_rate(bool use_delta_t)
{
    if (use_delta_t)
    {
    return get_delta_t() * lysis;}
    else
    {
        return lysis;
    }
}

double Automaton::get_lysogenization_rate(bool use_delta_t)
{
    if (use_delta_t)
    {
    return get_delta_t() * lysogenization;}
    else
    {
        return lysogenization;
    }
}

double Automaton::get_virion_degradation_rate()
{
    return get_delta_t() * virion_degradation;
}

double Automaton::get_ab_lysis_response(bool delta_t)
{
    if (delta_t)
    {
    return get_delta_t() * ab_lysis_response;}
    else
    {
        return ab_lysis_response;
    }
}

    //set data about bacteria
void Automaton::set_state(unsigned state_in)
  {
    state = state_in;
  }
void Automaton::set_bacterial_growth_rate(double growth_in)
{
    bacterial_growth_rate = growth_in;
}

void Automaton::set_bacterial_resistance(double b_resistance_in)
{
    bacterial_resistance = b_resistance_in; 
}

  //set data about plasmid
void Automaton::set_plasmid_conjugation_rate(double input)
{
    plasmid_conjugation_rate = input;
}

void Automaton::set_plasmid_loss_rate(double input)
{
    plasmid_loss_rate = input;
}

void Automaton::set_plasmid_resistance(double input)
{
    plasmid_resistance = input;
}

  //set data about phage
void Automaton::set_phage_resistance(double input)
{
    phage_resistance = input;
}

void Automaton::set_burst_size(unsigned input)
{
    burst_size = input;
}

void Automaton::set_infection_rate(double input)
{
    infection_rate = input;
}

void Automaton::set_lysis_rate(double input)
{
    lysis = input;
}

void Automaton::set_lysogenization_rate(double input)
{
    lysogenization = input;
}

void Automaton::set_virion_degradation_rate(double input)
{
    virion_degradation = input;
}

void Automaton::set_ab_lysis_response(double input)
{
    ab_lysis_response = input;
}