
from shutil import copyfile

# Open existing meta.yaml and another one
metayaml = open('meta.yaml')
outyaml = open('out.yaml', 'w')

# Find the build number, increment it, and write to the new yaml
found = False
for line in metayaml:
    if "number:" in line:
        found = True
        # For the line containing the build number, parse the number and increment
        elements = [e.strip() for e in line.split(":")]
        if not elements[1].isnumeric():
            raise ValueError("Build number is not parsable: {}".format(line))

        old_build_number = int(elements[1])
        new_build_number = old_build_number + 1

        # Write new build number to new yaml
        outyaml.write(line.replace(str(old_build_number), str(new_build_number)))
    else:
        # Write all other lines to new yaml
        outyaml.write(line)

if not found:
    raise Exception("Error incrementing the build number.")

# Clean up
metayaml.close()
outyaml.close()

# Replace original meta.yaml with the new one
copyfile('out.yaml', 'meta.yaml')
