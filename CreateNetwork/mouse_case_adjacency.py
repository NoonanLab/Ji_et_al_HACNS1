import sys, re

inp = sys.argv[1] if len(sys.argv) > 1 else "-"
fin = sys.stdin if inp == "-" else open(inp)

def mouse_case_token(tok: str) -> str:
    """Convert ALLCAPS gene-like tokens to Mouse style: first letter upper, rest lower (letters only).
       Preserve digits and separators; leave tokens with lowercase as-is."""
    tok = tok.strip()
    if not tok or tok == "NA":
        return tok
    # split on separators but keep them
    parts = re.split(r'([_\-/])', tok)
    out_parts = []
    for p in parts:
        if re.fullmatch(r'[_\-/]', p):
            out_parts.append(p)
            continue
        # If it's all caps/digits (typical human style), lower only letters after first letter
        if re.fullmatch(r'[A-Z0-9]+', p):
            m = re.match(r'([A-Z])([A-Z0-9]*)', p)
            if m:
                first, rest = m.group(1), m.group(2)
                rest = ''.join(ch.lower() if ch.isalpha() else ch for ch in rest)
                out_parts.append(first + rest)
            else:
                out_parts.append(p)
        else:
            # Mixed-case already (likely fine); still ensure first alpha is upper, rest alphas lower
            chars = list(p)
            first_alpha_idx = next((i for i,c in enumerate(chars) if c.isalpha()), None)
            if first_alpha_idx is not None:
                chars[first_alpha_idx] = chars[first_alpha_idx].upper()
                for j in range(first_alpha_idx+1, len(chars)):
                    if chars[j].isalpha():
                        chars[j] = chars[j].lower()
            out_parts.append(''.join(chars))
    return ''.join(out_parts)

header = fin.readline()
sys.stdout.write(header)  # "Source\tTargets\n"

for line in fin:
    if not line.strip():
        continue
    src, tgts = line.rstrip("\n").split("\t", 1)
    src2 = mouse_case_token(src)
    # split targets on commas, strip, transform, and join back with ", "
    tgt_list = [t.strip() for t in tgts.split(",")]
    tgt_list2 = [mouse_case_token(t) for t in tgt_list if t != ""]
    sys.stdout.write(f"{src2}\t{', '.join(tgt_list2)}\n")
