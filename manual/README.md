# The nwsbm manual

`nwsbm_manual.pdf` is the user manual and reference of nwsbm (title "Manual", no version number in the citation). The LaTeX source is `nwsbm_manual.tex` with one file per chapter (`ch*.tex`); the listings of the examples are `ex/out/*.txt`, the figures are PDF files in `fig/` made by `ex/ex05_plots.do`.

## Bumping the version (one-line change)

The manual's title and citation carry no version number; the title page states the series it describes (`\manualseries`, valid for 1.0.x). The exact version and date of the release that produced the listings are two more macros in `nwsbm_manual.tex`:

    \newcommand{\manualseries}{1.0}
    \newcommand{\pkgversion}{1.0.2}
    \newcommand{\pkgdate}{20 November 2026}

Change `\manualseries` only for a new minor or major release whose documented behaviour changed. For a patch release rebuild (`make -C manual pdf`) and rerun the examples if their output changed. (`tools/set_version.py` does not touch the manual; the sentence in `ch01_intro.tex` with the tag, `v1.0.0`, is changed by hand.)

## The DOI

`\manualdoi` is empty until the Zenodo record exists; an empty DOI prints nothing in the citation. `ZENODO.md` holds the metadata of the record and the steps to do.

## Building

    make -C manual pdf        # pdflatex, 3 passes
    make -C manual check      # undefined references, overfull boxes
    NWCORE=/path/to/nwcommands [NWROOT=/path/to/nwsbm] make -C manual examples   # reruns the do-files of ex/ in Stata 17

`NWROOT` defaults to this repository; set it to a checkout of a release tag to reproduce the listings of that release. The do-files `ex/ex01_karate.do` to `ex/ex06_compare.do` write one text log per listing into `ex/out` with `log using`; `ex/mklistings.py` removes the log header and footer lines. The data of the examples (Zachary, Lazega, Davis) are loaded from the nwcommands website (internet needed). The listings use fixed seeds; the timings in `ex/out/misc_timing.txt` depend on the computer and its load.
