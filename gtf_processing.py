import re
import sys

inp  = sys.argv[1]
out  = sys.argv[2]

with open(inp) as f, open(out, "w") as w:
    for line in f:
        if line.startswith("#"):
            w.write(line)
            continue
        parts = line.rstrip("\n").split("\t")
        attrs = parts[8]
        m = re.search(r'locus_tag\s+"([^"]+)"', attrs)
        if m:
            locus = m.group(1)
            attrs = re.sub(r'gene_id\s+"[^"]+"', f'gene_id "{locus}"', attrs)
        parts[8] = attrs
        w.write("\t".join(parts) + "\n")