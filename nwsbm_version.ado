*! version 1.1.0  23nov2026 nwsbm
* Copyright (C) 2026 Thomas U. Grund. Part of nwsbm; GNU GPL-3 (see LICENSE).
capture program drop nwsbm_version
program nwsbm_version, rclass
	version 14
	local __pkgver "1.1.0"
	local __pkgdate "23nov2026"
	di as text "nwsbm version `__pkgver', `__pkgdate'"
	di as text "  status                  : stochastic block models (Bernoulli, Poisson, Gaussian, degree-corrected, two-mode)"
	di as text "  license                 : GNU GPL-3"
	return local version "`__pkgver'"
	return local date "`__pkgdate'"
end
