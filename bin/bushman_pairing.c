#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <limits.h>

/* Members are unique ascending original key-row IDs, exactly as which(%in%).
 * Search the short side in a large hub instead of rehashing the hub per pair.
 * No candidate pair or supporting key row is discarded. */
static int common(const int *a, int na, const int *b, int nb, int *out) {
    int n = 0;
    if (!na || !nb) return 0;
    if (na > nb) {
        const int *p = a; a = b; b = p;
        int k = na; na = nb; nb = k;
    }
    if (na < nb / 8) {
        int lo = 0;
        for (int i = 0; i < na; ++i) {
            int left = lo, right = nb;
            while (left < right) {
                int mid = left + (right - left) / 2;
                if (b[mid] < a[i]) left = mid + 1; else right = mid;
            }
            lo = left;
            if (lo == nb) break;
            if (b[lo] == a[i]) { if (out) out[n] = a[i]; ++n; ++lo; }
        }
    } else {
        int i = 0, j = 0;
        while (i < na && j < nb) {
            if (a[i] < b[j]) ++i;
            else if (a[i] > b[j]) ++j;
            else { if (out) out[n] = a[i]; ++n; ++i; ++j; }
        }
    }
    return n;
}

static void check_members(SEXP values, SEXP ends) {
    if (TYPEOF(values) != INTSXP || TYPEOF(ends) != INTSXP ||
        XLENGTH(values) > INT_MAX || XLENGTH(ends) > INT_MAX)
        Rf_error("Compact members require integer vectors within IRanges limits");
    int previous = 0;
    for (R_xlen_t i = 0; i < XLENGTH(ends); ++i) {
        int end = INTEGER(ends)[i];
        if (end == NA_INTEGER || end < previous || end > XLENGTH(values))
            Rf_error("Invalid compact partition");
        for (int j = previous; j < end; ++j)
            if (INTEGER(values)[j] <= 0 ||
                (j > previous && INTEGER(values)[j] <= INTEGER(values)[j-1]))
                Rf_error("Compact members must be unique ascending positive row IDs");
        previous = end;
    }
    if (previous != XLENGTH(values)) Rf_error("Incomplete compact partition");
}

SEXP bushman_intersections(SEXP av, SEXP ae, SEXP bv, SEXP be, SEXP ai, SEXP bi) {
    check_members(av, ae); check_members(bv, be);
    if (TYPEOF(ai) != INTSXP || TYPEOF(bi) != INTSXP ||
        XLENGTH(ai) != XLENGTH(bi) || XLENGTH(ai) > INT_MAX)
        Rf_error("Invalid candidate locus slots");
    R_xlen_t pairs = XLENGTH(ai), total = 0;
    SEXP ends = PROTECT(Rf_allocVector(INTSXP, pairs));
    for (R_xlen_t i = 0; i < pairs; ++i) {
        int a = INTEGER(ai)[i], b = INTEGER(bi)[i];
        if (a == NA_INTEGER || b == NA_INTEGER || a < 1 || b < 1 ||
            a > XLENGTH(ae) || b > XLENGTH(be))
            Rf_error("Candidate locus slot out of bounds");
        int as = a == 1 ? 0 : INTEGER(ae)[a-2];
        int bs = b == 1 ? 0 : INTEGER(be)[b-2];
        total += common(INTEGER(av)+as, INTEGER(ae)[a-1]-as,
                        INTEGER(bv)+bs, INTEGER(be)[b-1]-bs, NULL);
        if (total > INT_MAX) Rf_error("Support exceeds IRanges integer partition limit; no data truncated");
        INTEGER(ends)[i] = (int)total;
        if ((i & 1023) == 0) R_CheckUserInterrupt();
    }
    SEXP values = PROTECT(Rf_allocVector(INTSXP, total));
    for (R_xlen_t i = 0; i < pairs; ++i) {
        int a = INTEGER(ai)[i], b = INTEGER(bi)[i];
        int as = a == 1 ? 0 : INTEGER(ae)[a-2];
        int bs = b == 1 ? 0 : INTEGER(be)[b-2];
        int out_start = i == 0 ? 0 : INTEGER(ends)[i-1];
        common(INTEGER(av)+as, INTEGER(ae)[a-1]-as,
               INTEGER(bv)+bs, INTEGER(be)[b-1]-bs, INTEGER(values)+out_start);
        if ((i & 1023) == 0) R_CheckUserInterrupt();
    }
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, values); SET_VECTOR_ELT(result, 1, ends);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("values"));
    SET_STRING_ELT(names, 1, Rf_mkChar("ends"));
    Rf_setAttrib(result, R_NamesSymbol, names);
    UNPROTECT(4);
    return result;
}

static const R_CallMethodDef calls[] = {
    {"bushman_intersections", (DL_FUNC) &bushman_intersections, 6},
    {NULL, NULL, 0}
};
void R_init_bushman_pairing(DllInfo *dll) {
    R_registerRoutines(dll, NULL, calls, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
