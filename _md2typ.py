"""Convert mathematics.md to Typst (.typ)."""
import re

PATH_MD = r"C:\Users\mccan\Code\avo_ibo\mathematics.md"
PATH_TYP = r"C:\Users\mccan\Code\avo_ibo\mathematics.typ"

with open(PATH_MD, encoding="utf-8") as f:
    text = f.read()

# ---------------------------------------------------------------------------
# Math conversion: LaTeX -> Typst, protect | inside math
# ---------------------------------------------------------------------------
def extract_braces(s, pos):
    """Return (content, end_pos) of balanced {...} starting at s[pos] (must be '{')."""
    depth = 1
    j = pos + 1
    while j < len(s):
        if s[j] == '{': depth += 1
        elif s[j] == '}':
            depth -= 1
            if depth == 0:
                return s[pos+1:j], j + 1
        j += 1
    return None, len(s)


def frac_replacer(s):
    r"""Replace \frac{a}{b} and \sqrt{a} (including nested braces)."""
    out = []
    i = 0
    while i < len(s):
        # Check for \sqrt{
        sqrt_m = re.match(r'\\sqrt\{', s[i:])
        if sqrt_m:
            i += sqrt_m.end()
            arg, i = extract_braces(s, i - 1)  # pass the { position
            if arg is not None:
                out.append(f'sqrt({arg})')
            continue

        # Check for \frac, \dfrac, \tfrac
        m = re.match(r'\\(?:d|t)?frac\{', s[i:])
        if m:
            i += m.end()  # after \frac{
            arg1, i = extract_braces(s, i - 1)  # i-1 is the { position
            if arg1 is None:
                continue
            # Second argument
            if i < len(s) and s[i] == '{':
                arg2, i = extract_braces(s, i)
                if arg2 is not None:
                    out.append(f'frac({arg1}, {arg2})')
            continue
        else:
            out.append(s[i])
            i += 1
    return ''.join(out)


def convert_math_source(src):
    """Convert LaTeX math body to Typst."""
    # 0. Vertical bars (BEFORE everything to avoid \lvert\tilde collision)
    src = re.sub(r'\\[lr]?vert(?!\w)', '|', src)
    src = re.sub(r'\\[lr]?Vert(?![a-zA-Z])', '||', src)

    # 1a. Merge \mathrm{word},more → \mathrm{word,more} so the comma is inside quotes
    for cmd in ['mathrm', 'operatorname']:
        while re.search(rf'\\{cmd}\{{([^}}]*)\}},(\w+)', src):
            src = re.sub(rf'\\{cmd}\{{([^}}]*)\}},(\w+)', rf'\\{cmd}{{\1,\2}}', src)

    # 1b. Font / function commands with braces
    src = re.sub(r'\\operatorname\{([^}]*)\}', r'"\1"', src)
    src = re.sub(r'\\mathrm\{([^}]*)\}', r'"\1"', src)
    src = re.sub(r'\\mathbf\{([^}]*)\}', r'bold(\1)', src)
    src = re.sub(r'\\mathcal\{([^}]*)\}', r'cal(\1)', src)
    src = re.sub(r'\\mathbb\{R\}', 'RR', src)
    src = re.sub(r'\\mathbb\{([^}]+)\}', r'mathbb(\1)', src)
    src = re.sub(r'\\tilde\{([^}]*)\}', r'tilde(\1)', src)
    src = re.sub(r'\\tilde\\([A-Za-z]+)', r'tilde(\1)', src)
    src = re.sub(r'\\hat\{([^}]*)\}', r'hat(\1)', src)
    src = re.sub(r'\\bar\{([^}]*)\}', r'bar(\1)', src)
    src = re.sub(r'\\text\{([^}]*)\}', r'"\1"', src)
    src = re.sub(r'\\mathtt\{([^}]*)\}', r'"\1"', src)
    src = re.sub(r'\\textbf\{([^}]*)\}', r'*\1*', src)
    src = re.sub(r'\\textit\{([^}]*)\}', r'_\1_', src)




    # \left / \right removal (BEFORE standalone cmds to avoid \le collision)
    src = re.sub(r'\\left', '', src)
    src = re.sub(r'\\right', '', src)

    # \frac, \dfrac, \tfrac (handles nested braces)
    src = frac_replacer(src)

    # 2. Standalone commands (longest-first for safety)
    # Greek → Unicode to avoid multi-letter identifier concatenation
    cmds = [
        (r'\varepsilon', 'epsilon'), (r'\varphi', 'varphi'),
        (r'\longrightarrow', '-->'), (r'\longleftarrow', '<--'),
        (r'\longmapsto', '|-->'), (r'\Longrightarrow', '=>'),
        (r'\Longleftarrow', '<='), (r'\leftrightarrow', '<->'),
        (r'\subseteq', 'subseteq'), (r'\supseteq', 'supseteq'),
        (r'\rightarrow', '->'), (r'\leftarrow', '<-'),
        (r'\Rightarrow', '=>'), (r'\Leftarrow', '<='),
        (r'\mapsto', '|->'), (r'\to', '->'), (r'\gets', '<-'),
        (r'\iff', '<=>'),
        (r'\partial', 'partial'), (r'\nabla', 'nabla'),
        (r'\infty', 'infty'), (r'\emptyset', 'emptyset'),
        (r'\forall', 'forall'), (r'\exists', 'exists'),
        (r'\approx', 'approx'), (r'\simeq', 'simeq'),
        (r'\sim', 'sim'), (r'\cong', 'cong'), (r'\equiv', 'equiv'),
        (r'\ne', '!='), (r'\neq', '!='),
        (r'\le', '<='), (r'\ge', '>='),
        (r'\lt', '<'), (r'\gt', '>'),
        (r'\ll', '<<'), (r'\gg', '>>'),
        (r'\subset', 'subset'), (r'\supset', 'supset'),
        (r'\notin', 'not in'), (r'\perp', 'perp'),
        (r'\parallel', 'parallel'),
        (r'\otimes', 'otimes'), (r'\oplus', 'oplus'),
        (r'\odot', 'odot'), (r'\circ', 'circ'),
        (r'\bullet', 'bullet'), (r'\cdot', 'cdot'),
        (r'\times', '\u00d7'), (r'\propto', 'propto'),
        (r'\alpha', 'alpha'), (r'\beta', 'beta'),
        (r'\gamma', 'gamma'), (r'\delta', 'delta'),
        (r'\theta', 'theta'), (r'\lambda', 'lambda'),
        (r'\mu', 'mu'), (r'\nu', 'nu'), (r'\xi', 'xi'),
        (r'\pi', 'pi'), (r'\rho', 'rho'),
        (r'\sigma', 'sigma'), (r'\tau', 'tau'),
        (r'\phi', 'phi'), (r'\chi', 'chi'),
        (r'\psi', 'psi'), (r'\omega', 'omega'),
        (r'\Gamma', 'Gamma'), (r'\Delta', 'Delta'),
        (r'\Theta', 'Theta'), (r'\Lambda', 'Lambda'),
        (r'\Pi', 'Pi'), (r'\Sigma', 'Sigma'),
        (r'\Phi', 'Phi'), (r'\Psi', 'Psi'), (r'\Omega', 'Omega'),
        (r'\sum', 'sum'), (r'\prod', 'prod'), (r'\int', 'int'),


        (r'\sin', 'sin '), (r'\cos', 'cos '), (r'\tan', 'tan '),
        (r'\arctan', 'arctan '), (r'\log', 'log '), (r'\exp', 'exp '),
        (r'\det', 'det '), (r'\arg', 'arg '),
        (r'\max', 'max '), (r'\min', 'min '),
        (r'\sup', 'sup '), (r'\inf', 'inf '), (r'\lim', 'lim '),
        (r'\colon', ':'), (r'\quad', 'quad'), (r'\qquad', 'quad quad'),
        (r'\,', ' '), (r'\:', ' '), (r'\;', ' '), (r'\!', ''),
        (r'\top', 'top'), (r'\bot', 'bot'),
        (r'\cdots', '..'), (r'\dots', '..'), (r'\ldots', '..'),
        (r'\in', 'in'),
    ]
    for latex, typst in sorted(cmds, key=lambda x: -len(x[0])):
        src = src.replace(latex, typst)

    # 3. Save literal \{ \} braces (set notation, etc.) with placeholders
    #    BEFORE step 4 strips backslashes, so we can distinguish from grouping braces
    src = src.replace(r'\{', '\x00LB\x00')
    src = src.replace(r'\}', '\x00RB\x00')

    # 3b. Bracket escapes (excluding braces which are handled above)
    src = src.replace(r'\[', '[')
    src = src.replace(r'\]', ']')
    src = src.replace(r'\langle', '\u27e8')
    src = src.replace(r'\rangle', '\u27e9')
    src = re.sub(r'\\big(g)?l(?!\w)', '', src)
    src = re.sub(r'\\big(g)?r(?!\w)', '', src)

    # 4. Strip remaining \cmd (literal-brace placeholders use \x00, not backslash → safe)
    src = re.sub(r'\\([A-Za-z]+)', r'\1', src)
    src = re.sub(r'\\tag\{[^}]*\}', '', src)

    # 5. Fix adjacent-letter concatenation in sub/superscript groups
    # e.g. _{mualpha} → _{mu alpha}, _{ij} → _{i j}
    GREEK_NAMES = ['alpha', 'beta', 'gamma', 'delta', 'epsilon',
        'varepsilon', 'zeta', 'eta', 'theta', 'iota', 'kappa', 'lambda',
        'mu', 'nu', 'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau',
        'upsilon', 'phi', 'varphi', 'chi', 'psi', 'omega',
        'Gamma', 'Delta', 'Theta', 'Lambda', 'Pi', 'Sigma', 'Phi',
        'Psi', 'Omega']

    # Multi-letter ASCII math keywords from cmds dict (excludes Greek names)
    MATH_KEYWORDS = sorted([
        'longrightarrow', 'longleftarrow', 'longmapsto', 'Longrightarrow',
        'Longleftarrow', 'leftrightarrow', 'subseteq', 'supseteq',
        'rightarrow', 'leftarrow', 'Rightarrow', 'Leftarrow',
        'mapsto', 'iff', 'partial', 'nabla', 'infty', 'emptyset',
        'forall', 'exists', 'approx', 'simeq', 'cong', 'equiv',
        'subset', 'supset', 'perp', 'parallel', 'otimes', 'oplus',
        'odot', 'circ', 'bullet', 'cdot', 'times', 'propto',
        'sin', 'cos', 'tan', 'arctan', 'log', 'exp', 'det', 'arg',
        'max', 'min', 'sup', 'inf', 'lim', 'in', 'notin',
        'colon', 'quad', 'top', 'bot',
        'cdots', 'dots', 'ldots',
    ], key=len, reverse=True)

    def fix_subsup_braces(m):
        content = m.group(2)
        # Protect quoted strings with non-alpha placeholders
        strings = []
        str_idx = [0]
        def save_str(m):
            idx = str_idx[0]
            str_idx[0] += 1
            strings.append(m.group(0))
            return f'\x00S{idx}\x00'
        content = re.sub(r'"[^"]*"', save_str, content)
        # Insert space before Greek names when preceded by a letter
        for name in sorted(GREEK_NAMES, key=len, reverse=True):
            content = re.sub(rf'(?<=[a-zA-Z])({re.escape(name)})(?![a-zA-Z])', r' \1', content)
        # Insert space after Greek names when followed by a letter (e.g. alphain → alpha in)
        for name in sorted(GREEK_NAMES, key=len, reverse=True):
            content = re.sub(rf'({re.escape(name)})(?=[a-zA-Z])', r'\1 ', content)
        # Protect Greek names with non-alpha placeholders
        for i, name in enumerate(GREEK_NAMES):
            content = content.replace(name, f'\x00{i}\x00')
        # Protect math keywords (multi-letter identifiers from LaTeX commands)
        for i, name in enumerate(MATH_KEYWORDS):
            content = content.replace(name, f'\x00K{i}\x00')
        # Insert space between adjacent Latin letters
        # (placeholders use \x00 + ... + \x00, no letters → not affected)
        content = re.sub(r'(?<=[a-zA-Z])(?=[a-zA-Z])', ' ', content)
        # Restore Greek names
        # Restore math keywords
        for i, name in enumerate(MATH_KEYWORDS):
            content = content.replace(f'\x00K{i}\x00', name)
        # Restore Greek names
        for i, name in enumerate(GREEK_NAMES):
            content = content.replace(f'\x00{i}\x00', name)
        # Restore strings
        for i, s in enumerate(strings):
            content = content.replace(f'\x00S{i}\x00', s)
        return m.group(1) + '(' + content + ')'

    src = re.sub(r'([_^])\{([^}]*)\}', fix_subsup_braces, src)

    # 5b. Strip remaining bare braces (LaTeX grouping; fix_subsup_braces already used ^( ) and _( ))
    src = src.replace('{', '').replace('}', '')

    # 5c. Restore literal braces from \{ \}
    src = src.replace('\x00LB\x00', '{')
    src = src.replace('\x00RB\x00', '}')

    # 6. Protect | (now from \vert restored or existing) with placeholder
    src = src.replace('|', '\x00PIPE\x00')

    return src


def replace_span(m):
    if m.group(2) is not None:       # $$...$$ → $...$ (single $ for display math in Typst 0.15)
        return '$\n' + convert_math_source(m.group(2)).strip() + '\n$'
    elif m.group(3) is not None:     # $...$
        return '$' + convert_math_source(m.group(3)) + '$'
    elif m.group(4) is not None:     # ```math ... ``` → display math
        return '$\n' + convert_math_source(m.group(4)).strip() + '\n$'
    return m.group(0)


# Combined pattern: $$..$$, $..$, ```math..``` ($$ before $)
pattern = re.compile(
    r'(\$\$(.*?)\$\$)'
    r'|(?<!\$)\$(?!\$)(.*?)(?<!\$)\$(?!\$)'
    r'|```math\n(.*?)```',
    re.DOTALL
)
text = pattern.sub(replace_span, text)
text = text.replace('```', '')

# ---------------------------------------------------------------------------
# Structural conversion
# ---------------------------------------------------------------------------
lines = text.split("\n")
output = []
i = 0
while i < len(lines):
    line = lines[i]

    h = re.match(r'^(#{1,4})\s+(.*)', line)
    if h:
        output.append("=" * (len(h.group(1)) - 1) + " " + h.group(2))
        i += 1
        continue

    if re.match(r'^-{3,}\s*$', line):
        output.append("#line()")
        i += 1
        continue

    if line.strip().startswith('|') and '|' in line:
        tbl = []
        while i < len(lines) and '|' in lines[i].strip():
            tbl.append(lines[i])
            i += 1
        data = [r for r in tbl if not re.match(r'^\|[-:\s|]+\|$', r)]
        if data:
            parsed = []
            for r in data:
                cells = [c.strip() for c in r.split('|')[1:-1]]
                parsed.append(cells)
            ncols = len(parsed[0])
            cs = ", ".join(["auto"] * ncols)
            output.append("#table(columns: (" + cs + "),")
            for ri, cells in enumerate(parsed):
                tcs = []
                for ci, cell in enumerate(cells):
                    tcs.append("[*" + cell + "*]" if ri == 0 else "[" + cell + "]")
                output.append("  " + ", ".join(tcs) + ",")
            output.append(")")
        continue

    output.append(line)
    i += 1

# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------
output = [l.replace('\x00PIPE\x00', '|') for l in output]
full = '\n'.join(output)
full = re.sub(r'\*\*(.+?)\*\*', r'*\1*', full, flags=re.DOTALL)
output = full.split('\n')

output = [re.sub(r'\[([^\]]+)\]\(([^)]+)\)', r'#link("\2")[\1]', l) for l in output]

preamble = """#set page(margin: 1in)
#set text(font: "New Computer Modern", size: 11pt)

"""

with open(PATH_TYP, "w", encoding="utf-8") as f:
    f.write(preamble)
    for l in output:
        f.write(l + "\n")

print(f"Written to {PATH_TYP}, {len(output)} lines")
