from glob import glob
import os
files = glob("*.pfout")
offset = 166 - 127
data = (("Liscenced Bed Capacity",134,2),
        ("Number of Closed Beds due to CP", 139, 2),
        ("Number of Closed Beds due to Staffing", 140, 2),
        ("Number of Bed moves due to Acuity", 141,2),
        ("Number of Mismatched Rooms MRSA Only", 149,2),
        ("Number of Mismatched Rooms VRE Only", 149,3),
        ("Number of Mismatched Rooms Both", 149,4),
        ("Number of Mismatched Rooms Total", 149,5),
        ("Number of Mismatched Rooms Due to Acuity", 144,2),
        ("Number of Mismatched Beds Due to Acuity", 145,2),
        ("Number of Mismatched Beds Due to Acuity (Low looking for High)", 146,2),
        ("Number of Mismatched Beds Due to Acuity (High looking for Low)", 147,2),
        )

for f in files:
    with open(f) as fread:
        lines = fread.readlines()

    hour = 0
    writefile = os.path.splitext(f)[0]+".extracted"
    with open(writefile, 'w') as fwrite:
        fwrite.write("time step")
        for output in data:
            label = output[0]
            fwrite.write("\t{}".format(label))
        
        while True:
            fwrite.write("\n")
            try:
                lines[128 + offset*hour]
                fwrite.write(str(hour))
            except IndexError:
                break
            for output in data:
                label, row, col = output
                fwrite.write("\t{}".format(lines[row+offset*hour].split("\t")[col].strip()))

            hour +=1
