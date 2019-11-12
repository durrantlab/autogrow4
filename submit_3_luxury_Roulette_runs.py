import os
import sys

#Submit Roulette luxury_big_run_mpi_Roulette.sh
sub_script_original = os.path.abspath("luxury_big_run_mpi_Roulette.sh")

output_dir = "/bgfs/jdurrant/jspiegel/Luxury_Run/Roulette_Lux/sub_scripts/"
if os.path.exists(output_dir)==False:
    os.mkdir(output_dir)

for i in [1,2,3]:
    printout = ""
    with open(sub_script_original, "r") as f:
        for line in f.readlines():
            if "#SBATCH --job-name=" in line:
                line = line.replace("_0\n","_{}\n".format(i))
            elif "#SBATCH --output=" in line:
                line = line.replace("0.conf.out\n","{}.conf.out\n".format(i))
            elif "for i in 1" in line:
                line = "for i in {}\n".format(i)
            printout = printout + line
    

    output_subscript = "{}luxury_big_run_mpi_Roulette_Run_{}.sh".format(output_dir,i)
    with open(output_subscript, "w") as f:
        f.write(printout) 

    command = "sbatch {}".format(output_subscript)
    print(command)
    os.system(command)

