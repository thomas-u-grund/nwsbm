#!/usr/bin/env python3
"""mklistings.py - the listings of the manual are the text logs ex/out/<name>.txt written by the do-files of ex/ (log using ..., text);
this script removes the header and footer lines that Stata adds, and the lines of the log commands, so that only the commands and their output remain."""
import glob, os, re
here = os.path.dirname(os.path.abspath(__file__))
drop = re.compile(r'^\s*(name:|log:|log type:|opened on:|closed on:)')
for f in sorted(glob.glob(os.path.join(here, 'out', '*.txt'))):
    lines = open(f, encoding='utf-8', errors='replace').read().split('\n')
    out = []
    skip_cont = False
    for l in lines:
        if skip_cont and l.startswith('> '):
            skip_cont = False
            continue
        skip_cont = False
        if drop.match(l):
            skip_cont = True
            continue
        if re.match(r'^\. (log close|log using)\b', l) or l.strip() == '.' and False: continue
        out.append(l.rstrip())
    # the rule that opens and closes a text log (as wide as the line size)
    def trim():
        while out and out[0] == '': out.pop(0)
        while out and out[-1] in ('', '.'): out.pop()
    trim()
    if out and len(out[0]) >= 80 and set(out[0]) == {'-'}: out.pop(0)
    if out and len(out[-1]) >= 80 and set(out[-1]) == {'-'}: out.pop()
    trim()
    open(f, 'w', encoding='utf-8').write('\n'.join(out) + '\n')
