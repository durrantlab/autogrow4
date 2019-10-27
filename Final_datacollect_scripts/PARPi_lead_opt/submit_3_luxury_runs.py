import os
import sys


sub_script_original = os.path.abspath("parpi_leadopt_MPI.sh")
output_dir = "/bgfs/jdurrant/jspiegel/PARPi_leadopt/sub_scripts/"

for i in [1,2,3]:
    printout = ""
    with open(sub_script_original,'r') as f:
        for line in f.readlines():
            if "#SBATCH --job-name=" in line:
                line = line.replace("_0\n","_{}\n".format(i))
            elif "#SBATCH --output=" in line:
                line = line.replace("0.conf.out\n","{}.conf.out\n".format(i))
            elif "for i in 1" in line:
                line = "for i in {}\n".format(i)
            printout = printout + line
    

    output_subscript = "{}PARPi_leadopt{}.sh".format(output_dir,i)
    with open(output_subscript, 'w') as f:
        f.write(printout) 

    command = "sbatch {}".format(output_subscript)
    print(command)
    os.system(command)

