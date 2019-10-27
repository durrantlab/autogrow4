
oldfile = "/home/jacob/Downloads/zinc15_available/pass_sanitize_all.smi_no_dup.smi"
new_file = "/home/jacob/Downloads/zinc15_available/pass_sanitize_all_no_dup2.smi"
counter = 0

line_counter= 0
with open(oldfile) as of:
    for lines in of.readlines():
        line_counter=line_counter+1




with open(oldfile) as of:
    with open(oldfile) as nf:
        while counter <line_counter:
            old_line = of.readlines()
            new_line = nf.readlines()
            counter =counter +1

            if old_line != new_line:
                print("")
                print("###############################")
                print("CHANGE IN line")
                print(old_line)
                print(new_line)
                print("###############################")

print(counter)
print("FINISHED!!!!!!!!")