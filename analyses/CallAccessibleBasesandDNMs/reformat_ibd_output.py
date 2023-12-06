#!/usr/bin/env python3
import sys
​
if len(sys.argv) != 2:
    print("Usage: ./reformat_ibd_output.py <input file>")
    sys.exit(1)
​
ibd_input_file = sys.argv[1]
​
numid2indiv = {
    "1079254": "P1",
    "1079255": "S1",
    "1079256": "C11",
    "1079257": "C12",
    "1079258": "P2",
    "1079259": "C21",
    "1079260": "C22",
    "1079261": "C23",
    "1079262": "P3",
    "1079263": "S3",
    "1079264": "C31",
    "1079265": "C32",
    "1079266": "P4",
    "1079267": "C41",
    "1079268": "C42",
}
output_lines = []
input_line_count = 0
with open(ibd_input_file, 'r') as f:
    for line in f:
        input_line_count += 1
        [s1, h1, s2, h2, chr_num, start, end, ibd_len] = line.rstrip().split("\t")
        # set haplotype to 'a' if h1 == 1 else 'b'
        s1_indiv, s2_indiv = numid2indiv[s1], numid2indiv[s2]
        sh1 = f"{s1_indiv}a" if h1 == "1" else f"{s1_indiv}b"
        sh2 = f"{s2_indiv}a" if h2 == "1" else f"{s2_indiv}b"
​
        output_lines.append([chr_num, start, end, sh1, sh2])
​
ibd_output_file = "reformat.ibd"
output_line_count = 0
with open(ibd_output_file, 'w') as f:
    for line in output_lines:
        output_line_count += 1
        f.write("\t".join(line) + "\n")
assert input_line_count == output_line_count, "Input and output line counts do not match"
