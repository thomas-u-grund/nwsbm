*! version 0.2.0  19nov2026 nwsbm
* Copyright (C) 2026 Thomas U. Grund. Part of nwsbm; GNU GPL-3 (see LICENSE).
capture program drop _nwsbm_libcheck
program _nwsbm_libcheck
	version 14
	capture mata: mata which SbmFit()
	if _rc {
		capture mata: mata mlib index
		capture mata: mata which SbmFit()
		if _rc {
			di as error "{p 0 0 2}nwsbm: the Mata library lnwsbm.mlib was not found. Reinstall the package with {bf:net install nwsbm, replace}.{p_end}"
			exit 3499
		}
	}
end
