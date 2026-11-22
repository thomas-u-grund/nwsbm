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
tempname fh
local out "dev/restarts/grid_part`part'.csv"
file open `fh' using "`out'", write replace
file write `fh' "dir,n,K,level,pin,pout,rep,Q_def,Q_fast,Q_ref,icl_def,icl_fast,icl_ref,j_def,j_fast,j_ref,ari_def_true,ari_fast_true,ari_ref_true,ari_def_fast,ari_def_ref,t_def,t_fast,t_ref" _n
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
        foreach v in def fast ref {
            local i = `i' + 1
            local opt ""
            if "`v'" == "fast" local opt "fast"
            if "`v'" == "ref" local opt "restarts(20)"
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
        foreach v in def fast ref {
            if `Q_`v'' < . {
                sbmt_ari m_`v' truth
                local at_`v' = r(ari)
            }
            else local at_`v' .
        }
        local adf .
        local adr .
        if `Q_def' < . & `Q_fast' < . {
            sbmt_ari m_def m_fast
            local adf = r(ari)
        }
        if `Q_def' < . & `Q_ref' < . {
            sbmt_ari m_def m_ref
            local adr = r(ari)
        }
        file write `fh' "`line',`Q_def',`Q_fast',`Q_ref',`icl_def',`icl_fast',`icl_ref',`j_def',`j_fast',`j_ref',`at_def',`at_fast',`at_ref',`adf',`adr',`t_def',`t_fast',`t_ref'" _n
    }
   }
  }
 }
}
file close `fh'
di as result "GRID DONE part `part'"
