import os

output_filename = "iladata.txt"

with open(output_filename, "w") as outfile:
    for filename in os.listdir("."):
        if filename.startswith("tmp_file") and os.path.isfile(filename):
            with open(filename, "r") as infile:
                lines = infile.readlines()[1:]
                outfile.writelines(lines)