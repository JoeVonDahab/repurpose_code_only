from pathlib import Path

fp = Path("myreceptor_targeted.gpf")

# Exact lines to drop entirely
drop_full_lines = {
    "map myreceptor_targeted.Si.map",
    "map myreceptor_targeted.CL.map",
    "map myreceptor_targeted.BR.map",
}

with fp.open("r", encoding="utf-8") as f:
    lines = f.readlines()

new_lines = []
for line in lines:
    stripped = line.strip()

    # Remove full map lines
    if stripped in drop_full_lines:
        continue

    # For ligand_types line, directly remove the " CL", " BR", " Si" substrings
    if line.startswith("ligand_types"):
        line = line.replace(" CL", "").replace(" BR", "").replace(" Si", "")

    new_lines.append(line)

with fp.open("w", encoding="utf-8") as f:
    f.writelines(new_lines)
