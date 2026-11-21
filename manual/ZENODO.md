# Zenodo record for the nwsbm manual

Metadata ready to paste into a **new** Zenodo record (the first for this manual; it is not a new version of the nwsaom or nwergm record). Nothing has been uploaded.

**DOI.** None yet. On zenodo.org start a new upload, press *Reserve DOI*, write the reserved all-versions DOI into `\manualdoi` in `nwsbm_manual.tex` (it is empty now, and an empty DOI prints nothing), rebuild with `make -C manual pdf` and upload the rebuilt PDF to the reserved record. The citation text always uses the all-versions DOI (master RELEASING.md B.10).

**Labelling.** The title and the citation of the manual contain no version number. The Zenodo version label is the series "X.Y" (here "1.0"; internal to Zenodo, not part of the citation), and the description says that the manual describes nwsbm X.Y.x. The PDF title page has one line saying which series it describes (`\manualseries`); the exact package version of the examples is `\pkgversion`. The software version is cited separately, as printed by `nwsbm version`.

**Patch releases.** The manual of the series 1.0 holds for 1.0.x; a patch release changes no documented option or number. Title and citation stay unchanged; no new Zenodo version is needed for a patch. A new Zenodo version only for 1.1 or 2.0 with changed documented behaviour (then change `\manualseries`).

## Metadata

**Title.** nwsbm: Stochastic block models for networks in Stata. Manual

**Version.** 1.0 (the series label, internal to Zenodo: the manual describes nwsbm 1.0.x; not part of the citation)

**Publication date.** The day of the upload.

**Resource type.** Publication, sub-type Software documentation (as the nwsaom and nwergm records).

**Creators.** Grund, Thomas U. (RWTH Aachen University), ORCID 0000-0001-6751-5429 (as on the other records; remove it in the form if the author does not want it).

**License.** Creative Commons Attribution 4.0 International (CC-BY-4.0), as the other manual records. The code stays under GPL-3.0-only.

**Access.** Open.

**Description** (plain text, one paragraph):

> This is the user manual and reference of nwsbm, series 1.0 (examples produced with version 1.0.2 and the nwcommands core 3.33.0), a Stata package for stochastic block models of networks: the actors are sorted into blocks, the distribution of a tie depends only on the blocks of the two actors, and the number of blocks is chosen by the integrated classification likelihood (ICL). nwsbm builds on the nwcommands core package and follows the R packages blockmodels and sbm (variational EM, objective, ICL, model names) and greed (exact ICL, degree-corrected models, greedy search). It fits Bernoulli, Poisson and Gaussian models for undirected and directed one-mode networks with missing ties, the degree-corrected Poisson model, the latent block model for two-mode networks, and dyadic covariates with standard errors of the coefficients. The manual describes the installation, the models and the estimator, the syntax and every option, the selection of the number of blocks, the exact ICL and the greedy search, covariates, two-mode networks, the postestimation commands and stored results, the plots (blocked matrix, expected ties, memberships, block graph, ICL path), the comparison of models, the validation against sbm, blockmodels and greed as recorded in the test suite of the package, and the limits. Worked examples use the karate club, the law firm of Lazega and the Southern women of Davis, Gardner and Gardner. nwsbm requires the nwcommands core package (version 3.29.0 or later) and Stata 17 or later and is meant for researchers who work in Stata. It is not a replacement for the R packages, which remain the reference implementations, and it is not affiliated with or endorsed by their authors. The software is free software under the GNU General Public License, version 3. The package and the manual were prepared with the help of an AI coding assistant working under the author's direction, as stated in the manual. Please cite this manual (all-versions DOI) and report the exact version printed by nwsbm version.

**Keywords.** Stata; nwcommands; nwsbm; stochastic block model; latent block model; degree-corrected block model; social network analysis; integrated classification likelihood; variational EM; blockmodels; sbm; greed

**Language.** English

**Related identifiers.** (relation, identifier, type)

| Relation | Identifier | Type |
|---|---|---|
| Is supplement to | https://github.com/thomas-u-grund/nwsbm/releases/tag/v1.0.2 | URL |
| Is supplement to | https://github.com/thomas-u-grund/nwsbm | URL |
| Requires | https://github.com/thomas-u-grund/nwcommands | URL |
| Is related to (manual of the sibling package nwsaom) | 10.5281/zenodo.23623807 | DOI |
| Is related to (manual of the sibling package nwergm) | 10.5281/zenodo.23661331 | DOI |

Confirm before the upload that the GitHub release v1.0.2 is published (if the package has a later 1.0.x release by then, use its tag). The two sibling links are optional.

**Other fields.** Funding: none. Communities: none.

## Steps for the owner

1. Read the PDF once (`manual/nwsbm_manual.pdf`); the examples were produced with Stata 17 on macOS with version 1.0.2 of the package (tag v1.0.2) and core 3.33.0.
2. On zenodo.org create the new record, reserve the DOI, put it into `\manualdoi`, rebuild (`make -C manual pdf`), upload the PDF and fill the fields above.
3. Commit the rebuilt PDF and the figures (`git add -f manual/nwsbm_manual.pdf manual/fig/*.pdf`).
4. Afterwards: the DOI goes into the citation texts (README of the package, `CITATION.cff`, the website page *How to cite* of nwsbm and the manual link of its menu); none of these is changed yet.
