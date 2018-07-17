from __future__ import division
from pathos.helpers import mp
import numpy as np

def weighted_high_median(a, wts):
    N = len(a)
    wtotal = 0
    wdiscardedlow = 0

    for i in xrange(N):
        wtotal += wts[i]

    nn = N
    while True:
        assert (nn > 0 and len(a) == nn)

        trial = sorted(a)[int(nn/2)]
        # Count up the weight to the left of and at the trial point.
        # Weight to the right of it isn't needed
        wleft = wtrial = 0
        for i in xrange(nn):
            if a[i] < trial:
                wleft += wts[i]
            elif a[i] == trial:
                wtrial += wts[i]

        if 2*(wdiscardedlow + wleft) > wtotal:
            # Trial value is too high
            ncandidates = 0
            #for i = 1:nn
            for i in xrange(nn):
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
            for i in xrange(nn):
                if a[i] > trial:
                    a[ncandidates] = a[i]
                    wts[ncandidates] = wts[i]
                    ncandidates += 1
            nn = ncandidates
            wdiscardedlow += wleft+wtrial

        a=a[:nn]
        wts=wts[:nn]

def qn(data):
    # sort data
    data = np.sort(data)

    n = len(data)
    h = int(n/2) + 1
    k = int(h*(h-1)/2)

    left = np.arange(n+1,1,-1)
    right = np.full(n,n, dtype= np.int64)

    work = np.zeros(n) # dtype = np.float64
    weight = np.zeros(n, np.int64)
    P = np.zeros(n, np.int64)
    Q = np.zeros(n, np.int64)

    jhelp = int((n*(n+1))/2)
    knew = k+jhelp
    nL = jhelp
    nR = n*n
    found = False
    Qn = 0*data[0]

    while (nR-nL) > n:
        j = 1
        for i in xrange(1,n,1):
            if left[i] <= right[i]:
                weight[j-1] = right[i] - left[i] + 1
                jhelp = left[i] + int(weight[j-1]/2)
                work[j-1] =data[i] - data[n+1-jhelp-1]
                j += 1

        trial = weighted_high_median(work[:j-1], weight[:j-1])

        j=0
        for i in xrange(n-1, -1, -1):
            while (j < n) and (data[i]-data[n-j-1] < trial):
                j += 1
            P[i] = j

        j = n+1
        for i in xrange(n):
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
            found = True
            break

    if found == False:
        j=1
        for i in xrange(1,n,1):
            if left[i] <= right[i]:
                for jj in xrange(left[i], right[i]+1, 1):
                    work[j-1] = data[i]-data[n-jj]
                    j += 1

        Qn = sorted(work[:j])[knew-nL-1]

    if n<10:
        nscale = [0, .399, .994, .512, .844, .611, .857, .669, .872][n-1]
    elif n%2 == 1:
        nscale = n/(n+1.4)
    else:
        nscale = n/(n+3.8)

    Qn = Qn*2.2219*nscale

    return Qn

def multiple_testing_correction(pval_dict):
    import statsmodels.api as sm

    nodepair_list = pval_dict.keys()
    qval_list = sm.stats.multipletests([pval_dict[nodepair] for nodepair in nodepair_list], method='fdr_bh')[1].tolist()
    qval_dict = {}
    for _, (i,j) in enumerate(nodepair_list):
        qval_dict[(j,i)] = qval_dict[(i,j)] = qval_list[_]
    return qval_dict

def get_cluster_size_distribution(clusterid_to_taxa):
    clusterlen_to_frequency = {}
    for clusterid, taxa in clusterid_to_taxa.items():
        try:
            clusterlen_to_frequency[len(taxa)] += 1
        except:
            clusterlen_to_frequency[len(taxa)] = 1

    return [i for j in [[clusterlen] * frequency for clusterlen, frequency in clusterlen_to_frequency.items()] for i in j]

def summary_stats(clusterid_to_taxa, master_leafpair_to_distance, master_nodepair_to_dist, clusterlen_distribution, statsfname, treefname, total_clustered_count, total_taxa_count, min_cluster_size, fdr_cutoff, gamma, hytest_method, dispersion_method, qval_determination, within_cluster_limit, solution_index, clean_up_status, pc_input):

    #def summary_stats(clusterid_to_taxa, master_leafpair_to_distance, master_nodepair_to_dist, clusterlen_distribution, statsfname, treefname, total_clustered_count, total_taxa_count, min_cluster_size, fdr_cutoff, gamma, hytest_method, dispersion_method, qval_determination, within_cluster_limit, solution_index, pc_input):

    """
    Summary stats of within- and inter-clade divergence
    """
    import itertools

    cluster_to_mean_and_med_dist = {}

    for clusterid, taxa in clusterid_to_taxa.items():
        pairwise_leaf_distance_distribution = [master_leafpair_to_distance[(leaf_x, leaf_y)] for leaf_x, leaf_y in itertools.combinations(taxa, 2)]
        cluster_to_mean_and_med_dist[clusterid] = (np.mean(pairwise_leaf_distance_distribution), np.median(pairwise_leaf_distance_distribution))

    # within cluster statistics
    # mean pwdist of all clusters
    mean_dist = [_[0] for _ in cluster_to_mean_and_med_dist.values()]

    # median pwdist of all clusters
    median_dist = [_[-1] for _ in cluster_to_mean_and_med_dist.values()]

    # separation between clusters
    intercluster_dist = [master_nodepair_to_dist[i][j] for i,j in itertools.combinations(clusterid_to_taxa.keys(), 2)]
    median_intercluster = np.median(intercluster_dist)
    mad_intercluster = np.median([abs(x-median_intercluster) for x in intercluster_dist])*1.4826

    # cluster size distributions
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
                   mu_icd, sd_intercluster_dist, median_intercluster, mad_intercluster, min(intercluster_dist), max(intercluster_dist)]

        output.write('{}\t{}\t{}\t{}\t'
                     '{}\t{}\t{}\t{}\t{}\t{}\t'
                     .format(treefname, min_cluster_size, fdr_cutoff, gamma,
                             hytest_method, dispersion_method, qval_determination, within_cluster_limit, solution_index, clean_up_status))

        #output.write('{}\t{}\t{}\t{}\t'
        #             '{}\t{}\t{}\t{}\t{}\t'
        #             .format(treefname, min_cluster_size, fdr_cutoff, gamma,
        #                     hytest_method, dispersion_method, qval_determination, within_cluster_limit, solution_index))

        if pc_input == False:
            output.write('na\t')
        else:
            prior_taxa = [i for j in pc_input.values() for i in j]
            clustered_taxa = [i for j in clusterid_to_taxa.values() for i in j]
            output.write('{}\t'.format(len(set(prior_taxa)&set(clustered_taxa))/len(prior_taxa)))

        output.write('{}\r\n'.format('\t'.join(map(str, results))))

    return (coverage, mu_pwd, mu_icd)