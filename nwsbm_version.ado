*! version 0.3.0  19nov2026 nwsbm
* Copyright (C) 2026 Thomas U. Grund. Part of nwsbm; GNU GPL-3 (see LICENSE).
capture program drop nwsbm_version
program nwsbm_version, rclass
	version 14
	local v "0.1.0"
	local d "19nov2026"
	di as text "nwsbm version `v', `d'"
	di as text "  status                  : development (Bernoulli stochastic block model)"
	di as text "  license                 : GNU GPL-3"
	return local version "`v'"
	return local date "`d'"
end
