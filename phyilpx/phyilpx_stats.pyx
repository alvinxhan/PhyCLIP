from __future__ import division
import numpy as np
import cython
import itertools
cimport numpy as np

IF UNAME_SYSNAME == "Windows":
    ctypedef long long _int64
ELSE:
    ctypedef long _int64

cdef extern from "math.h":
    float exp "exp" (float)
    float sqrt "sqrt" (float)
    float fabs "fabs" (float)
    _int64 floor "floor" (float)

def hypotest(np.ndarray data1, np.ndarray data2, _int64 method):
    if method == 1:
        return kuiper(data1, data2)
    else:
        return ks_2samp(data1, data2)

@cython.boundscheck(False)
cdef float _kuiper_dist(float q):

    cdef float EPS1 = 0.001
    cdef float EPS2 = 0.00000001

    cdef float a2 = -2.*q*q
    cdef float term = 0.
    cdef float sum = 0.
    cdef float termbf = 0.
    cdef float dist = 1.

    cdef _int64 j

    for j in range(1, 100001):
        term = 2*((4*j*j*q*q)-1)*exp(a2*j*j)
        sum += term

        if fabs(term) <= EPS1*termbf or fabs(term) <= EPS2*sum:
            dist = sum
            break

        termbf = fabs(term)
    return dist

@cython.boundscheck(False)
cdef float kuiper(float [:] data1, float [:] data2):
    cdef _int64 j1 = 0
    cdef _int64 j2 = 0
    cdef _int64 n1 = len(data1)  # number of elements in data1
    cdef _int64 n2 = len(data2)  # number of elements in data2
    cdef float eff_n = sqrt((n1*n2)/(n1+n2))
    cdef float fn1 = 0.
    cdef float fn2 = 0.
    cdef float d1, d2, dt, lam, d = 0.

    while j1 < n1 and j2 < n2:
        d1 = data1[j1]
        d2 = data2[j2]

        if d1 <= d2:
            j1 += 1
            fn1 = j1/n1

        if d2 <= d1:
            j2 += 1
            fn2 = j2/n2

        dt = fabs(fn2-fn1)
        if dt > d:
            d = dt

    lam = (eff_n + 0.155 + 0.24/eff_n)*d


    return _kuiper_dist(lam)

@cython.boundscheck(False)
cdef float ks_2samp(float [:] data1, float [:] data2):
    from scipy.stats import kstwobign

    cdef _int64 n1 = len(data1) # number of elements in data1
    cdef _int64 n2 = len(data2) # number of elements in data2

    cdef np.ndarray data_all = np.concatenate((data1, data2)) # concatenate both data arrays and sort - this is basically the array of all possible values

    cdef np.ndarray cdf1 = np.searchsorted(data1,data_all,side='right')/n1
    cdef np.ndarray cdf2 = np.searchsorted(data2,data_all,side='right')/n2

    cdef float d = np.amax(np.absolute(cdf1-cdf2))
    # Note: d absolute not signed distance
    cdef float en = sqrt((n1*n2)/(n1+n2))

    cdef float prob
    try:
        prob = kstwobign.sf((en + 0.12 + 0.11 / en) * d)
    except:
        prob = 1.0
    return prob


@cython.boundscheck(False)
cdef float weighted_high_median(np.ndarray a, np.ndarray wts):
    cdef _int64 N = len(a)
    cdef _int64 wtotal = 0
    cdef _int64 wdiscardedlow = 0

    cdef _int64 i

    for i in range(N):
        wtotal += wts[i]

    cdef _int64 nn = N
    cdef float trial
    cdef _int64 wleft
    cdef _int64 wtrial
    cdef _int64 ncandidates

    while True:
        assert (nn > 0 and len(a) == nn)

        #trial = np.sort(a)[int(nn/2)]
        trial = np.sort(a)[floor(nn/2)]
        # Count up the weight to the left of and at the trial point.
        # Weight to the right of it isn't needed
        wleft = wtrial = 0
        for i in range(nn):
            if a[i] < trial:
                wleft += wts[i]
            elif a[i] == trial:
                wtrial += wts[i]

        if 2*(wdiscardedlow + wleft) > wtotal:
            # Trial value is too high
            ncandidates = 0
            #for i = 1:nn
            for i in range(nn):
                if a[i] < trial:
                    a[ncandidates] = a[i]
                    wts[ncandidates] = wts[i]
                    ncandidates += 1
            nn = ncandidates

        elif 2*(wdiscardedlow + wleft + wtrial) > wtotal:
            # Trial value is just right
            return trial

        else:
            # Trial value is too low
            ncandidates = 0
            #for i = 1:nn
            for i in range(nn):
                if a[i] > trial:
                    a[ncandidates] = a[i]
                    wts[ncandidates] = wts[i]
                    ncandidates += 1
            nn = ncandidates
            wdiscardedlow += wleft + wtrial

        a=a[:nn]
        wts=wts[:nn]

def qn(data):
    return _qn(np.array(data, dtype = np.float32))

@cython.boundscheck(False)
cdef float _qn(np.ndarray data):

    # sort data
    data = np.sort(data)

    cdef _int64 n = len(data)
    #cdef _int64 h = int(n/2) + 1
    cdef _int64 h = floor(n/2) + 1
    #cdef _int64 k = int(h*(h-1)/2)
    cdef _int64 k = floor(h*(h-1)/2)

    cdef np.ndarray left = np.arange(n+1,1,-1)
    cdef np.ndarray right = np.full(n, n, dtype = np.int64)

    cdef np.ndarray work = np.zeros(n, dtype = np.float32) # dtype = np.float32
    cdef np.ndarray weight = np.zeros(n, dtype = np.int64)
    cdef np.ndarray P = np.zeros(n, dtype = np.int64)
    cdef np.ndarray Q = np.zeros(n, dtype = np.int64)

    #cdef _int64 jhelp = int((n*(n+1))/2)
    cdef _int64 jhelp = floor((n*(n+1))/2)
    cdef _int64 knew = k+jhelp
    cdef _int64 nL = jhelp
    cdef _int64 nR = n*n
    cdef _int64 found = 0
    cdef float Qn = 0*data[0]

    cdef _int64 j
    cdef _int64 i
    cdef _int64 sumP
    cdef _int64 sumQ
    cdef float trial
    cdef float nscale

    while (nR-nL) > n:
        j = 1
        for i in range(1,n,1):
            if left[i] <= right[i]:
                weight[j-1] = right[i] - left[i] + 1
                #jhelp = left[i] + int(weight[j-1]/2)
                jhelp = left[i] + floor(weight[j-1]/2)
                work[j-1] =data[i] - data[n+1-jhelp-1]
                j += 1

        trial = weighted_high_median(work[:j-1], weight[:j-1])

        j=0
        for i in range(n-1, -1, -1):
            while (j < n) and (data[i]-data[n-j-1] < trial):
                j += 1
            P[i] = j

        j = n+1
        for i in range(n):
            while data[i]-data[n-j+2-1] > trial: # 55
                j -= 1
            Q[i] = j

        sumP = sum(P)
        sumQ = sum(Q)-n # 60

        if knew <= sumP:
            right[:] = P[:]
            nR = sumP
        elif knew > sumQ:
            left[:] = Q[:]
            nL = sumQ
        else:
            Qn = trial
            found = 1
            break

    if found == 0:
        j=1
        for i in range(1,n,1):
            if left[i] <= right[i]:
                for jj in range(left[i], right[i]+1, 1):
                    work[j-1] = data[i]-data[n-jj]
                    j += 1

        Qn = np.sort(work[:j])[knew-nL-1]

    if n<10:
        nscale = [0, .399, .994, .512, .844, .611, .857, .669, .872][n-1]
    elif n%2 == 1:
        nscale = n/(n+1.4)
    else:
        nscale = n/(n+3.8)

    return Qn*2.2219*nscale


def summary_stats(object clusterid_to_taxa, np.ndarray master_nodepair_to_dist, object clusterlen_distribution, str statsfname, str treefname, _int64 total_clustered_count, _int64 total_taxa_count, _int64 min_cluster_size, float fdr_cutoff, float gamma, str hytest_method, str dispersion_method, str qval_determination, int force, float within_cluster_limit, _int64 solution_index, str sol_ver, str sol_mode, object pc_input):

    """
    Summary stats of within- and inter-clade divergence
    """
    cdef object cluster_to_mean_and_med_dist = {}
    cdef object taxa, pairwise_leaf_distance_distribution, mean_dist, median_dist, intercluster_dist

    cdef _int64 clusterid

    for clusterid, taxa in clusterid_to_taxa.items():
        pairwise_leaf_distance_distribution = [master_nodepair_to_dist[(leaf_x, leaf_y)] for leaf_x, leaf_y in itertools.combinations(taxa, 2)]
        cluster_to_mean_and_med_dist[clusterid] = (np.mean(pairwise_leaf_distance_distribution), np.median(pairwise_leaf_distance_distribution))

    # within cluster statistics
    # mean pwdist of all clusters
    cdef float median_intercluster, mad_intercluster
    mean_dist = [_[0] for _ in cluster_to_mean_and_med_dist.values()]

    # median pwdist of all clusters
    median_dist = [_[-1] for _ in cluster_to_mean_and_med_dist.values()]

    # separation between clusters
    intercluster_dist = [master_nodepair_to_dist[(i, j)] for i,j in itertools.combinations(clusterid_to_taxa.keys(), 2)]
    median_intercluster = np.median(intercluster_dist)
    mad_intercluster = np.median([abs(x-median_intercluster) for x in intercluster_dist])*1.4826

    # cluster size distributions
    cdef float median_cluster_size, mad_cluster_size
    median_cluster_size = np.median(clusterlen_distribution)
    mad_cluster_size = np.median([abs(x-median_cluster_size) for x in clusterlen_distribution])*1.4826

    with open(statsfname, 'a') as output:
        try:
            sd_mean_dist = np.std(mean_dist, ddof=1)
        except:
            sd_mean_dist = 'nan'

        try:
            sd_median_dist = np.std(median_dist, ddof=1)
        except:
            sd_median_dist = 'nan'

        try:
            sd_intercluster_dist = np.std(intercluster_dist, ddof=1)
        except:
            sd_intercluster_dist = 'nan'

        try:
            sd_cluster_size = np.std(clusterlen_distribution, ddof=1)
        except:
            sd_cluster_size = 'nan'

        coverage = total_clustered_count/total_taxa_count
        mu_pwd = np.mean(mean_dist)
        mu_icd = np.mean(intercluster_dist)

        results = [total_clustered_count, total_taxa_count, coverage, len(clusterid_to_taxa),
                   np.mean(clusterlen_distribution), sd_cluster_size, median_cluster_size, mad_cluster_size, min(clusterlen_distribution), max(clusterlen_distribution),
                   mu_pwd, sd_mean_dist, np.mean(median_dist), sd_median_dist, min(mean_dist), max(mean_dist),
                   mu_icd, sd_intercluster_dist, median_intercluster, mad_intercluster, min(intercluster_dist), max(intercluster_dist), sol_ver]

        """
        results = [total_clustered_count, total_taxa_count, coverage, len(clusterid_to_taxa),
                   np.mean(clusterlen_distribution), sd_cluster_size, median_cluster_size, mad_cluster_size, min(clusterlen_distribution), max(clusterlen_distribution),
                   mu_pwd, sd_mean_dist, np.mean(median_dist), sd_median_dist, min(mean_dist), max(mean_dist),
                   mu_icd, sd_intercluster_dist, median_intercluster, mad_intercluster, min(intercluster_dist), max(intercluster_dist), sol_ver, sol_mode]
        """

        output.write('{}\t{}\t{:.3f}\t{:2f}\t'
                     '{}\t{}\t{}\t{}\t{}\t{}\t'
                     .format(treefname, min_cluster_size, fdr_cutoff, gamma,
                             hytest_method, dispersion_method, qval_determination, force, within_cluster_limit, solution_index))

        if pc_input == False:
            output.write('na\t')
        else:
            prior_taxa = [i for j in pc_input.values() for i in j]
            clustered_taxa = [i for j in clusterid_to_taxa.values() for i in j]
            output.write('{}\t'.format(len(set(prior_taxa)&set(clustered_taxa))/len(prior_taxa)))

        output.write('{}\r\n'.format('\t'.join(map(str, results))))

    return (coverage, mu_pwd, mu_icd)