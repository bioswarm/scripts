# Build mapping to new B-chain residue numbers
mapping = {}
new_resseq = 1
for key in residues_LM:
    mapping[key] = new_resseq
    new_resseq += 1

new_lines = []

for line in lines:
    record = line[:6].strip()
    if record in ("ATOM", "HETATM"):
        if len(line) < 27:
            new_lines.append(line)
            continue
        chain = line[21]
        try:
            resseq = int(line[22:26])
        except ValueError:
            new_lines.append(line)
            continue
        icode = line[26]
        key = (chain, resseq, icode)
        if chain in ("L", "M") and key in mapping:
            new_chain = "B"
            new_res = mapping[key]
            line = (
                line[:21]
                + new_chain
                + f"{new_res:4d}"
                + line[26:]
            )
        elif chain == "R":
            line = line[:21] + "A" + line[22:]
    elif record == "TER":
        # Update chain ID in TER lines if long enough
        if len(line) >= 22:
            chain = line[21]
            if chain in ("L", "M"):
                line = line[:21] + "B" + line[22:]
            elif chain == "R":
                line = line[:21] + "A" + line[22:]
    new_lines.append(line)

out_path = "/mnt/data/Toll-Spz-new_merged.pdb"
with open(out_path, "w") as f:
    f.writelines(new_lines)

out_path, os.path.getsize(out_path)
