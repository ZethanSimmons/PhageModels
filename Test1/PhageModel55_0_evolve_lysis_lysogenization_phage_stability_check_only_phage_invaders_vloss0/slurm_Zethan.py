import sys
import os
import os.path
import stat
import shutil
import subprocess
import time

"""
MEMO: This experiment is to check what lysogenization and lysis rates phages spontaneously evolve. It also runs a new experiment at the end to check whether starting with the evolved rate at t0 makes any difference. It has influx from only phage-bearing genomine configurations
"""


home = os.getenv("HOME")

needed_files = [home + '/Thesis/PhageModel45_0_evolve_lysis_lysogenization_phage_stability_check_only_phage_invaders_CORRECTED/demo', home + '/Thesis/PhageModel45_0_evolve_lysis_lysogenization_phage_stability_check_only_phage_invaders_CORRECTED/source.tar.gz', home + '/Thesis/PhageModel45_0_evolve_lysis_lysogenization_phage_stability_check_only_phage_invaders_CORRECTED/slurm_Zethan.py']

def main():
    base_dir = '/nesi/nobackup/uoa04131/data/PhageModel/PhageModel55_0_evolve_lysis_lysogenization_phage_stability_check_only_phage_invaders_vloss0/'

    command = """#!/bin/bash -e

#SBATCH --job-name=550v0_{pops:}  # job name (shows up in the queue)
#SBATCH --account=uoa04131      # Project Account
#SBATCH --time=167:00:00         # Walltime (HH:MM:SS)
#SBATCH --mem=1900      # memory/cpu (in MB)
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --cpus-per-task=1       # number of cores per task (e.g. OpenMP) 36 cores/node
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --output=out-%j      # %x and %j are replaced by job name and ID
#SBATCH --error=err-%j

srun demo --population {pops:} --s_antibiotic_time_periodicity {ab_period:} --s_antibiotic_time_limit {ab_lim:} --s_the_seed {seed:} --a_antibiotics {antibiotics:} --a_v_resistance {vres:} --a_p_resistance {pres:} --a_b_resistance {bres:} --a_lysogenization_rate {lysg:} --a_lysogenization_rate_001 {lysg1:} --a_lysis_rate {lys:} --a_lysis_rate_001 {lys1:} --a_conjugation {conj:} --a_conjugation_010 {conj1:} --s_random_mixing {rand_mix:} --s_death_lysis_switch {dlys:} --a_b_death {deth:} --a_v_diffusion {v_diffus:} --s_run1 {run1:} --s_run2 {run2:} --a_v_loss_rate 0.0

"""
    #Note that all influx rates have already been set in main.cpp, due to the more convenient interface there
    populations = (901, 101, 109)
    ab_treatments = ((1,10), (1, 100), (1,1000))
    ab_strengths = (0.3,)
    conjugation_rates = (0.0,)
    lysis_rates = (0.005,)
    lysogenization_rates = (0.85,)
    random_mixing_switches = (0,1) 
    death_rates = (0.12,) #0.02, 0.06
    death_lysis_switches = (0,1)
    v_diffusion_rates = (0.18,)
    label_count = 0

    for population in populations:
        for ab_treatment in ab_treatments:
            for ab_strength in ab_strengths:
                for conj_rate0 in conjugation_rates:
                    for lys_rate0 in lysis_rates:
                        for lysg_rate0 in lysogenization_rates:
                            for rmixswitch in random_mixing_switches:
                                for death_rate in death_rates:
                                    for death_lysis_switch in death_lysis_switches:
                                        for v_diffusion_rate in v_diffusion_rates:
                                            label_count+=1
                                            per_experiment_dir = base_dir + "pop" + str(population) + "_AB_periodicity" + str(ab_treatment[1]) + "_lys" + str(lys_rate0) + "_lysg" + str(lysg_rate0) + "_conj" + str(conj_rate0) + "_rmix" + str(rmixswitch) + "_deth" + str(death_rate) + "_dlys" + str(death_lysis_switch) + "_vdiff" + str(v_diffusion_rate)
                                            runsim(per_experiment_dir, command.format(pops=population, ab_period=ab_treatment[1], ab_lim = ab_treatment[0], seed=uqnum(), antibiotics = ab_strength, vres = ab_strength, pres=ab_strength, bres = ab_strength, lysg = lysg_rate0, lysg1 = lysg_rate0, lys = lys_rate0, lys1 = lys_rate0, conj=conj_rate0, conj1 = conj_rate0, rand_mix=rmixswitch, dlys=death_lysis_switch, deth=death_rate, v_diffus=v_diffusion_rate, run1=0, run2 = 1, count=label_count))





def static_var(varname, value):
    def decorate(func):
        setattr(func, varname, value)
        return func
    return decorate

@static_var("this_is_first_call", True)
def uqnum():
    """This function returns a unique integer as a string upon every call. The number can be used as a random seed."""
    if uqnum.this_is_first_call:
        uqnum.num = int(time.time()*10**5) % 2**31 # This will not work 'forever'.
        uqnum.this_is_first_call = False
    else:
        uqnum.num += 1

    return '{here:d}'.format(here=uqnum.num)


def runsim(directory, command):
    # Check if a directory where a simulation is run already exists. If
    # the directory does not exist, make it
    if not os.path.exists(directory):
        try:
            subprocess.run(['/bin/mkdir', '-p', directory], check=True)
        except subprocess.CalledProcessError as err:
            print('ERROR:', err, file=sys.stderr)
            exit(1)
    # else:
    #     sys.stderr.write('ERROR:' + directory + ' already exists\n')
    #     exit(1)

    # Copy necessary files to the directory
    for file in needed_files:
        dest_file = os.path.join(directory, os.path.basename(file))
        if os.path.exists(dest_file):
            backup_file = dest_file + ".bak"
            shutil.copy2(dest_file, backup_file)
        try:
            shutil.copy2(file, directory)
        except shutil.Error as err:
            print('ERROR:', err, file=sys.stderr)
            exit(1)

    # Write command.sh
    with open(os.path.join(directory, 'command.sh'), 'w') as fout:
        fout.write(command)

    exec_mode = stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR
    os.chmod(os.path.join(directory, 'command.sh'), exec_mode)

    # Run command.sh
    try:
        result = subprocess.run(['sbatch', 'command.sh'], check=True, cwd=directory, stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as err:
        print('ERROR:', err, file=sys.stderr)
        exit(1)
    else:
        with open(os.path.join(directory, 'job_id'), 'w') as fout:
            print(result.stdout.decode(sys.stdout.encoding).strip(),file=fout)

if __name__ == '__main__':
    main()
