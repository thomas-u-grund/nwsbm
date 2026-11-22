* Validation grid for the default number of restarts: default (10 extra random starts) vs fast vs restarts(20).
* Run from the repository root: NWSBM_PART=<0..3> NWSBM_NPARTS=4 [NWSBM_REPS=10] stata -e do dev/restarts/grid.do ; output dev/restarts/grid_part<k>.csv
do cscripts/nwsbm_testlib.do
mata:
mata set matastrict off
void gsim(real scalar n, real scalar K, real scalar pin, real scalar pout, real scalar dir)
{
	real colvector z
	real matrix P, A, U
	z = mod((1::n), K) :+ 1
	z = z[order(uniform(n, 1), 1)]
	P = (z * J(1, n, 1) :== J(n, 1, 1) * z') * (pin - pout) :+ pout
	U = uniform(n, n)
	A = (U :< P)
	if (!dir) A = (lowertriangle(A, 0) + lowertriangle(A, 0)')
	_diag(A, 0)
	st_matrix("__g_A", A)
	st_matrix("__g_z", z)
}
end
local part : env NWSBM_PART
local nparts : env NWSBM_NPARTS
local reps : env NWSBM_REPS
if "`part'" == "" local part 0
if "`nparts'" == "" local nparts 1
if "`reps'" == "" local reps 10
local only500 : env NWSBM_N500
tempname fh
local out "dev/restarts/grid2`only500'_part`part'.csv"
file open `fh' using "`out'", write replace
file write `fh' "dir,n,K,level,pin,pout,rep,Q_r0,Q_r2,Q_r4,icl_r0,icl_r2,icl_r4,j_r0,j_r2,j_r4,ari_r0_true,ari_r2_true,ari_r4_true,ari_r0_r4,ari_r2_r4,t_r0,t_r2,t_r4" _n
local cell 0
foreach dir in 0 1 {
 foreach n in 50 100 200 500 {
  foreach K in 2 3 5 8 {
   foreach lev in easy mid hard linux {
    local cell = `cell' + 1
    if mod(`cell' + floor(`cell' / 4) + floor(`cell' / 16), `nparts') != `part' continue
    if "`lev'" == "easy" { 
        local pin 0.40
        local pout 0.05 
    }
    if "`lev'" == "mid" {
        local pin 0.25
        local pout 0.10
    }
    if "`lev'" == "hard" {
        local pin 0.16
        local pout 0.10
    }
    if "`lev'" == "linux" {
        local pin 0.25
        local pout 0.08
    }
    local mb = min(`K' + 3, 10)
    if ("`only500'" == "" & `n' == 500) | ("`only500'" != "" & `n' != 500) continue
    local nrep = cond(`n' == 500, min(`reps', 3), `reps')
    forvalues rep = 1/`nrep' {
        set seed `=100000 * `dir' + 1000 * `n' + 100 * `K' + `rep' + 7 * strlen("`lev'")'
        mata: gsim(`n', `K', `pin', `pout', `dir')
        clear
        quietly set obs `n'
        forvalues j = 1/`n' {
            quietly gen double __a`j' = .
        }
        mata: st_store(., ., st_matrix("__g_A"))
        quietly gen truth = .
        mata: st_store(., "truth", st_matrix("__g_z"))
        if `dir' quietly nwset __a*, name(g) directed replace
        else quietly nwset __a*, name(g) undirected replace
        local line "`dir',`n',`K',`lev',`pin',`pout',`rep'"
        local i 0
        foreach v in r0 r2 r4 {
            local i = `i' + 1
            local opt ""
            if "`v'" == "r0" local opt "restarts(0)"
            if "`v'" == "r2" local opt "restarts(2)"
            if "`v'" == "r4" local opt "restarts(4)"
            timer clear 1
            timer on 1
            capture quietly nwsbm g, maxblocks(`mb') seed(`rep') threads(1) `opt' generate(m_`v') replace
            local rc = _rc
            timer off 1
            quietly timer list 1
            local t_`v' = r(t1)
            if `rc' {
                local Q_`v' .
                local icl_`v' .
                local j_`v' .
            }
            else {
                local Q_`v' = e(Q)
                local icl_`v' = e(ICL)
                local j_`v' = e(J)
            }
        }
        foreach v in r0 r2 r4 {
            if `Q_`v'' < . {
                sbmt_ari m_`v' truth
                local at_`v' = r(ari)
            }
            else local at_`v' .
        }
        local a02 .
        local a24 .
        if `Q_r0' < . & `Q_r4' < . {
            sbmt_ari m_r0 m_r4
            local a02 = r(ari)
        }
        if `Q_r2' < . & `Q_r4' < . {
            sbmt_ari m_r2 m_r4
            local a24 = r(ari)
        }
        file write `fh' "`line',`Q_r0',`Q_r2',`Q_r4',`icl_r0',`icl_r2',`icl_r4',`j_r0',`j_r2',`j_r4',`at_r0',`at_r2',`at_r4',`a02',`a24',`t_r0',`t_r2',`t_r4'" _n
    }
   }
  }
 }
}
file close `fh'
di as result "GRID DONE part `part'"
