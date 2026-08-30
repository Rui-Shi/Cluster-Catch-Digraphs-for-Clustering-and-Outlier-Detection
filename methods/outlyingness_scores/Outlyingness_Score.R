# Calculate the vicinity density of each points
# data: the entire data set or cluster
# the radius of each point
# d: dimensionality
Vic_Den = function(data, R, d){
  n = dim(data)[1]
  ddatax = as.matrix(dist(data))
  if(is.null(n)) {n=1; ddatax=matrix(0,nrow=1,ncol=1)} # the case when there is only one observation
  Den = sapply(1:n, function(x){
    size = length(which(ddatax[x,]<=R[x]))
    # return(size/(R[x]^d))
    # return((size/R[x]^d)^(1/d))  # overflows: R^d = Inf for R > exp(709/d)
    #                              # (~13.3 at d=274), giving Den = 0 -> OOS NaN,
    #                              # IOS Inf. See revision_experiments/tr2/
    #                              # test_outlyingness_density.R (2026-07-16).
    return(size^(1/d)/R[x])       # algebraically identical, numerically stable
  })
  return(Den)
}

# Calculate the Outbound Outlying-ness Score (OOS) of each points
# data: the entire data set or cluster
# the radius of each point
# d: dimensionality
OOS = function(data, R, d){
  n = dim(data)[1]
  ddatax = as.matrix(dist(data))
  if(is.null(n)) {n=1; ddatax=matrix(0,nrow=1,ncol=1)} # the case when there is only one observation
  diag(ddatax) = Inf
  Den = Vic_Den(data, R, d)
  scores = sapply(1:n,function(x){
    out_nei.index = which(ddatax[x,]<=R[x])
    score = mean(Den[out_nei.index])/Den[x]
  })
  return(scores)
}

# Calculate the Inbound Outlying-ness Score (OOS) of each points
# data: the entire data set or cluster
# the radius of each point
# d: dimensionality
IOS = function(data, R, d){
  n = dim(data)[1]
  ddatax = as.matrix(dist(data))
  if(is.null(n)) {n=1; ddatax=matrix(0,nrow=1,ncol=1)} # the case when there is only one observation
  Den = Vic_Den(data, R, d)
  scores = sapply(1:n,function(x){
    in_nei.index = which(ddatax[x,]<=R)
    score = 1/sum(Den[in_nei.index])
  })
  return(scores)
}

# standardization using MADN, with the manuscript's degenerate-case fallback
# chain MADN -> SD -> 0 (CCDwScores.tex line 625, item B3): when MADN = 0 the
# denominator is replaced by the sample SD (median centering unchanged); if SD
# also vanishes (or is undefined, i.e. a singleton cluster), all scores are 0.
# See revision_experiments/tr2/test_std_madn.R.
# x: a set of sample point
std_MADN = function(x){
  if(mad(x)!=0){
    # s = abs((x-median(x))/mad(x))
    s = (x-median(x))/mad(x) # allow negative scores
  } else {
    sdx = sd(x)
    if(!is.na(sdx) && sdx!=0){
      s = (x-median(x))/sdx  # MADN=0 fallback: replace denominator with SD
    } else {
      s = rep(0,length(x))   # MADN=SD=0 (or singleton): degenerate cluster
    }
  }
  return(s)
}

# ---------------------------------------------------------------------------
# Tie-breaking for outlyingness scores (manuscript Eq. (9)-(10), Sec. 2.2)
# ---------------------------------------------------------------------------
# Manuscript setting: within a cluster C, sort the IOS values ascending and
# suppose a block of m points is tied,
#
#   IOS_(1) <= ... <= IOS_(k-1) <= IOS_(k) = ... = IOS_(k+m-1) <= IOS_(k+m) <= ... <= IOS_(n_c),
#
# bracketed below by the next distinct smaller value IOS_(k-1) and above by
# the next distinct larger value IOS_(k+m). Eq. (10) linearizes the tied
# block using each point's vicinity density rho_{x_i}, i = 1..m:
#
#   IOS~(x_i) := IOS_(k+m) - (IOS_(k+m) - IOS_(k-1)) * rho_{x_i} / sum_{j=1}^{m} rho_{x_j}.
#
# Direction property (manuscript prose, Sec. 2.2): the adjusted scores lie
# STRICTLY between IOS_(k-1) and IOS_(k+m), and within the tied block a
# SMALLER vicinity density yields a LARGER adjusted score (rho_{x_i} -> 0
# maps to the upper bracket IOS_(k+m); rho_{x_i} -> sum(rho) maps to the
# lower bracket IOS_(k-1)).
#
# 'legacy' exists only to reproduce previously published numbers: the
# original inline loop (present verbatim, four times, in NNCCD_OOS/IOS and
# RKCCD_OOS/IOS) computes
#
#   IOS_(k-1) + (IOS_(k+m) - IOS_(k-1)) * rho_{x_i} / sum(rho),
#
# the mirror image of Eq. (10) about the midpoint of the bracket -- it
# reverses the within-block order the manuscript specifies. break_ties()
# preserves that behaviour bit-exactly under mode = "legacy" so old results
# can still be regenerated; new work should use the corrected mode = "eq10".
#
# Known divergences between the legacy loop and the manuscript, and their
# disposition here:
#   (a) legacy uses 'break' instead of 'next' when a guard fails, so the
#       first unusable tie block aborts ALL remaining tie-breaking for that
#       call. FIXED in 'eq10' (uses 'next' via the outer while-loop).
#   (b) legacy matches tie blocks by relative tolerance
#       abs(v - x) / x < 1e-4 rather than exact equality. FIXED in 'eq10'
#       (exact equality; ties here arise from identical arithmetic).
#   (c) legacy runs tie-breaking over the whole input vector, while the
#       manuscript defines it within a cluster C. OUT OF SCOPE for this
#       patch -- 'eq10' keeps the same whole-vector scope as legacy; the
#       calling function is responsible for passing per-cluster vectors
#       where the manuscript requires that scope.
#   (d) legacy is applied to OOS as well as IOS, while the manuscript
#       defines tie-breaking only for IOS. OUT OF SCOPE for this patch --
#       'eq10' does not restrict which score type it is called on; that
#       remains a call-site decision.
#   (e) legacy does not guard non-finite bracket values (OOS can be +Inf).
#       FIXED in 'eq10' (a block is skipped, left tied, whenever its lower
#       or upper bracket is non-finite, or its bracket is degenerate, or
#       its density sum is non-finite/non-positive).
#
# scores: numeric vector of scores (any order; ties broken in place, values
#         returned in the same order as the input).
# vd:     numeric vector of vicinity densities, same length/order as scores.
# mode:   "eq10" (correct, default) | "legacy" (bit-exact old behaviour) |
#         "none" (no tie-breaking). Call sites take the default from
#         getOption("ccd.tiebreak", "eq10") rather than hard-coding it, so a
#         whole run can be switched between modes via options().
break_ties = function(scores, vd, mode = c("eq10", "legacy", "none")){
  mode = match.arg(mode)
  if(mode=="none"){
    return(scores)
  } else if(mode=="legacy"){
    return(break_ties_legacy(scores, vd))
  } else {
    return(break_ties_eq10(scores, vd))
  }
}

# Pre-existing behaviour, reproduced bit-exactly (this is the mirror-image
# formula described above, with divergences (a), (b), (e) also intact).
# Kept only so previously published numbers can be regenerated; do not use
# for new work.
break_ties_legacy = function(scores, vd){
  score_vd = cbind(scores,vd)
  order_index = order(scores)
  score_vd = score_vd[order_index,]
  frequency = table(scores)
  repeated_values = as.numeric(names(frequency[frequency > 1]))
  for(x in repeated_values){
    if(x==0){
      index = which(score_vd[,1]==x)
    } else {
      index = which(abs(score_vd[,1]-x)/x<0.0001)
    }
    if(length(index)<2) break
    index_max = max(index)
    index_min = min(index)
    if(index_min==1){
      s = score_vd[,1][1]
      e = score_vd[,1][index_max+1]
      diff = e - s
    } else if(index_max==length(scores)){
      s = score_vd[,1][index_min-1]
      e = score_vd[,1][index_max]
      diff = e - s
    } else {
      s = score_vd[,1][index_min-1]
      e = score_vd[,1][index_max+1]
      diff = e - s
    }
    if(is.na(diff)) break
    den_sum = sum(score_vd[,2][index])
    score_vd[,1][index] = sapply(index, function(x){
      new_s = s + diff*score_vd[,2][x]/den_sum
      return(new_s)
    })
    scores = score_vd[order(order_index),1]
  }
  return(scores)
}

# Correct implementation of manuscript Eq. (10). Fixes divergences (a), (b),
# (e) above; (c) and (d) are call-site scope decisions and are unchanged by
# this patch (see the comment block above).
break_ties_eq10 = function(scores, vd){
  n = length(scores)
  if(n<2) return(scores)
  order_index = order(scores)
  s = scores[order_index]
  v = vd[order_index]

  i = 1
  while(i <= n){
    # extend the run of exactly-equal values starting at i; NA/NaN never
    # equal anything (including themselves) under ==, so isTRUE(NA) is
    # FALSE and a block can never extend across or start on an NA/NaN --
    # such points are left untouched as length-1 "blocks" below.
    j = i
    while(j < n && isTRUE(s[j+1]==s[i])) j = j+1
    if(j > i){
      lo = if(i>1) s[i-1] else s[i]   # IOS_(k-1); own tied value if block is at the bottom
      hi = if(j<n) s[j+1] else s[j]   # IOS_(k+m); own tied value if block is at the top
      vd_block = v[i:j]
      den_sum = sum(vd_block)
      if(is.finite(lo) && is.finite(hi) && hi>lo &&
         is.finite(den_sum) && den_sum>0){
        s[i:j] = hi - (hi-lo)*vd_block/den_sum
      }
      # else: leave this block tied (non-finite bracket, degenerate bracket
      # e.g. the whole vector tied, or non-finite/non-positive density
      # sum) and move on -- this is the 'next' semantics called for above.
    }
    i = j+1
  }
  scores[order_index] = s
  return(scores)
}
