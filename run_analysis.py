from __future__ import annotations
import argparse, json, math, os, gc
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from collatz_compression import (
    accelerated_step, collatz_raw_metrics, iid_anticompression_probability,
    js_divergence, mutual_information_discrete, geometric_cramer_rate_lambda,
    residual_kernel_extrema, residual_bounds_fixed_total, near_extremal_count_bound, extremal_valuation_words, valuation_word_residue, valuation_word_density_among_odds, count_word_in_odd_progression, valuation_word_nondescend_threshold, valuation_word_witness_count, valuation_word_witness_ratio, witness_word_count_upper_bound, strict_supercritical_witness_exponent, exact_witness_word_count_small, LN2, LN3
)

SEED = 20260819


def ensure_dirs(root: Path):
    for p in [root/'results', root/'figures', root/'tables', root/'data']:
        p.mkdir(parents=True, exist_ok=True)


def odd_starts(base:int, n:int):
    s = base if base % 2 else base + 1
    return [s + 2*i for i in range(n)]


def wilson_interval(k:int, n:int, z:float=1.959963984540054):
    if n == 0: return (math.nan, math.nan)
    ph=k/n; den=1+z*z/n
    ctr=(ph+z*z/(2*n))/den
    half=z*math.sqrt((ph*(1-ph)+z*z/(4*n))/n)/den
    return max(0,ctr-half), min(1,ctr+half)


def _spectrum_chunk(args):
    starts,max_m=args
    out=np.empty((len(starts),max_m),dtype=np.int16)
    for idx,s in enumerate(starts):
        x=int(s)
        for j in range(max_m):
            x,k=accelerated_step(x)
            out[idx,j]=k
    return out

def spectrum_analysis(nstarts:int, ms:list[int], max_m:int, root:Path, base:int=10**30+1, workers:int=1):
    starts=odd_starts(base,nstarts)
    if workers>1:
        chunks=[c.tolist() for c in np.array_split(np.asarray(starts,dtype=object),workers) if len(c)]
        with Pool(processes=workers) as pool:
            parts=pool.map(_spectrum_chunk,[(c,max_m) for c in chunks])
        Kword=np.vstack(parts)
    else:
        Kword=_spectrum_chunk((starts,max_m))
    Kcum=np.cumsum(Kword,dtype=np.int32,axis=1)

    rows=[]
    for m in ms:
        C=Kcum[:,m-1]*LN2-m*LN3
        lam=C/m
        hits=int(np.sum(C<0)); anti=hits/nstarts
        lo,hi=wilson_interval(hits,nstarts)
        piid=iid_anticompression_probability(m)
        rows.append({'m':m,'mean_lambda':float(lam.mean()),'sd_lambda':float(lam.std(ddof=1)),
                     'anti_count':hits,'anti_prob_actual':anti,'anti_ci95_low':lo,'anti_ci95_high':hi,
                     'anti_prob_iid':piid,'actual_to_iid_ratio':anti/piid if piid>0 else math.nan,
                     'tail_rate_actual':-math.log(anti)/m if anti>0 else math.inf,
                     'tail_rate_iid':-math.log(piid)/m if piid>0 else math.inf,
                     'mean_K_over_m':float(Kcum[:,m-1].mean()/m)})
    sdf=pd.DataFrame(rows); sdf.to_csv(root/'tables'/'compression_summary.csv',index=False)

    cap=10; kvals=Kword.ravel()
    actual=np.array([(kvals==k).mean() for k in range(1,cap)] + [(kvals>=cap).mean()])
    baseline=np.array([2.0**(-k) for k in range(1,cap)] + [2.0**(-(cap-1))]); baseline/=baseline.sum()
    js=js_divergence(actual,baseline)
    pd.DataFrame({'k':[str(k) for k in range(1,cap)]+[f'>={cap}'],'actual':actual,'iid_geometric':baseline}).to_csv(root/'tables'/'valuation_distribution.csv',index=False)

    Kc=np.minimum(Kword,cap); dep=[]
    for lag in range(1,9):
        a=Kc[:,:max_m-lag].ravel(); b=Kc[:,lag:].ravel()
        dep.append({'lag':lag,'pearson_corr':float(np.corrcoef(a,b)[0,1]),'mutual_information_nats':mutual_information_discrete(a,b)})
    depdf=pd.DataFrame(dep); depdf.to_csv(root/'tables'/'dependence.csv',index=False)

    # Exact rare-event comparison.
    plt.figure(figsize=(7,4.5))
    plt.plot(sdf.m,sdf.anti_prob_actual,marker='o',label='Collatz ensemble')
    plt.plot(sdf.m,sdf.anti_prob_iid,marker='s',label='Exact residue-density law')
    plt.yscale('log'); plt.xlabel('Block length m (odd steps)'); plt.ylabel(r'$P(C_m<0)$'); plt.legend(); plt.tight_layout()
    plt.savefig(root/'figures'/'anti_compression_probability.pdf'); plt.close()

    # Tail rate at lambda=0.
    plt.figure(figsize=(7,4.5))
    plt.plot(sdf.m,sdf.tail_rate_actual,marker='o',label='Empirical')
    plt.plot(sdf.m,sdf.tail_rate_iid,marker='s',label='Exact density tail')
    plt.xlabel('Block length m'); plt.ylabel(r'$-m^{-1}\log P(\Lambda_m<0)$'); plt.legend(); plt.tight_layout()
    plt.savefig(root/'figures'/'tail_rate_zero.pdf'); plt.close()

    # Empirical finite-m rate curves, left of the mean, against Cramer rate.
    rate_rows=[]; rate_ms=[m for m in [16,32,64,128] if m<=max_m]
    grid=np.linspace(-0.12,0.24,55)
    for m in rate_ms:
        lam=(Kcum[:,m-1]*LN2-m*LN3)/m
        for z in grid:
            p=float(np.mean(lam<=z))
            if p>=1/nstarts and p<0.45:
                rate_rows.append({'m':m,'lambda':z,'empirical_rate':-math.log(p)/m,'cdf':p,'cramer_rate':geometric_cramer_rate_lambda(float(z))})
    rdf=pd.DataFrame(rate_rows); rdf.to_csv(root/'tables'/'rate_function.csv',index=False)
    plt.figure(figsize=(7,4.8))
    for m,g in rdf.groupby('m'):
        plt.plot(g['lambda'],g['empirical_rate'],marker='.',linewidth=1,label=f'm={m}')
    xg=np.linspace(-0.12,0.24,250); yg=[geometric_cramer_rate_lambda(float(x)) for x in xg]
    plt.plot(xg,yg,linewidth=2,label='Exact-density Cramer rate')
    plt.xlabel(r'$\lambda$'); plt.ylabel(r'$-m^{-1}\log P(\Lambda_m\leq\lambda)$'); plt.legend(ncol=2); plt.tight_layout()
    plt.savefig(root/'figures'/'rate_function_collapse.pdf'); plt.close()

    # Distribution at the longest scale.
    m=max(ms); lam=(Kcum[:,m-1]*LN2-m*LN3)/m
    plt.figure(figsize=(7,4.5)); plt.hist(lam,bins=80,density=True)
    plt.axvline(math.log(4/3),linestyle='--',label=r'Exact-density mean $\log(4/3)$')
    plt.xlabel(r'$\Lambda_m$'); plt.ylabel('Density'); plt.legend(); plt.tight_layout(); plt.savefig(root/'figures'/'compression_spectrum.pdf'); plt.close()

    plt.figure(figsize=(7,4.5)); xx=np.arange(len(actual)); w=.38
    plt.bar(xx-w/2,actual,w,label='Collatz'); plt.bar(xx+w/2,baseline,w,label='Geometric(1/2)')
    plt.xticks(xx,[str(k) for k in range(1,cap)]+[f'>={cap}']); plt.xlabel(r'$k=\nu_2(3n+1)$'); plt.ylabel('Probability'); plt.legend(); plt.tight_layout()
    plt.savefig(root/'figures'/'valuation_distribution.pdf'); plt.close()

    plt.figure(figsize=(7,4.5)); plt.plot(depdf.lag,depdf.mutual_information_nats,marker='o'); plt.xlabel('Lag'); plt.ylabel('Mutual information (nats)'); plt.tight_layout(); plt.savefig(root/'figures'/'dependence_mi.pdf'); plt.close()

    # Realized ordering effects at m=32 (or longest <=32 if necessary).
    mo=32 if max_m>=32 else max_m
    pref=np.cumsum(Kword[:,:mo],dtype=np.int32,axis=1)
    # S=1 + sum_{i=1}^{m-1} 2^{K_i}/3^i, stable via exp(C_i)
    Cs=pref[:,:mo-1]*LN2-np.arange(1,mo)*LN3
    log_terms=np.clip(Cs,-700,700)
    kernels=1.0+np.exp(log_terms).sum(axis=1)
    Ktot=pref[:,mo-1]
    targets=sorted(pd.Series(Ktot).value_counts().head(6).index.tolist())
    orows=[]
    for K in targets:
        mask=Ktot==K; lo,hi=residual_kernel_extrema(mo,int(K))
        orows.append({'m':mo,'K':int(K),'count':int(mask.sum()),'realized_kernel_min':float(kernels[mask].min()),'realized_kernel_median':float(np.median(kernels[mask])),
                      'realized_kernel_max':float(kernels[mask].max()),'sharp_all_words_min':lo,'sharp_all_words_max':hi,'max_fraction_of_sharp':float(kernels[mask].max()/hi)})
    odf=pd.DataFrame(orows); odf.to_csv(root/'tables'/'exceptional_word_ordering.csv',index=False)

    np.savez_compressed(root/'data'/'spectrum_valuations.npz', valuations=Kword, odd_interval_start=np.array([starts[0]],dtype=object), odd_interval_end=np.array([starts[-1]],dtype=object))
    meta={'n_odd_starts':nstarts,'odd_interval_start':starts[0],'odd_interval_end':starts[-1],'max_m':max_m,'js_divergence_k_vs_iid':js,'seed':SEED}
    (root/'results'/'spectrum_meta.json').write_text(json.dumps(meta,indent=2))
    # Release the large word matrices before trajectory analyses begin.
    del Kword, Kcum, Kc, pref, Cs, kernels, Ktot
    gc.collect()
    return sdf,depdf,rdf,odf,js


def trajectory_metrics_range(start:int, ncount:int, max_steps:int=10000):
    rows=[]
    for n in range(start,start+ncount):
        st,mx,ok=collatz_raw_metrics(n,max_steps)
        x=n
        while x%2==0 and x>1: x//=2
        n0=x; c=0.; d=0.; odd_steps=0; max_res=0.; end_c=0.; end_r=0.
        while x!=1 and odd_steps<2000:
            x,k=accelerated_step(x); odd_steps+=1; c += k*LN2-LN3; d=max(d,-c)
            rr=math.log(x/n0)+c; max_res=max(max_res,rr); end_c=c; end_r=rr
        rows.append((n,st,mx,ok,odd_steps,d,max_res,end_c,end_r,math.log(mx/n),math.log(n)))
    return pd.DataFrame(rows,columns=['n','total_stopping_time','max_excursion','reached_1','odd_steps','compression_deficit','max_residual','terminal_compression','terminal_residual','log_excursion_ratio','log_n'])


def fit_ols(train:pd.DataFrame,test:pd.DataFrame,features:list[str],target='log_excursion_ratio'):
    means=train[features].mean(); stds=train[features].std(ddof=0).replace(0,1)
    Xtr=(train[features]-means)/stds; Xte=(test[features]-means)/stds
    A=np.column_stack([np.ones(len(Xtr)),Xtr.to_numpy()]); y=train[target].to_numpy()
    beta=np.linalg.lstsq(A,y,rcond=None)[0]
    pred_tr=A@beta; pred_te=np.column_stack([np.ones(len(Xte)),Xte.to_numpy()])@beta
    def r2(yy,pp): return 1-float(np.sum((yy-pp)**2))/float(np.sum((yy-yy.mean())**2))
    return beta,r2(y,pred_tr),r2(test[target].to_numpy(),pred_te),means,stds,pred_te


def trajectory_analysis(nraw:int,root:Path,max_steps=10000):
    df=trajectory_metrics_range(1,nraw,max_steps); df.to_csv(root/'data'/'trajectory_metrics.csv.gz',index=False,compression='gzip')
    summary=pd.DataFrame([
        ['N',nraw],['reached_1',int(df.reached_1.sum())],['max_stopping_time',int(df.total_stopping_time.max())],
        ['argmax_stopping_time',int(df.loc[df.total_stopping_time.idxmax(),'n'])],['median_stopping_time',float(df.total_stopping_time.median())],
        ['max_excursion',int(df.max_excursion.max())],['argmax_excursion',int(df.loc[df.max_excursion.idxmax(),'n'])],['median_excursion',float(df.max_excursion.median())]
    ],columns=['metric','value']); summary.to_csv(root/'tables'/'trajectory_summary.csv',index=False)

    rho,p=spearmanr(df.compression_deficit,df.log_excursion_ratio)
    rng=np.random.default_rng(SEED); idx=rng.permutation(len(df)); cut=int(.7*len(df)); tr=df.iloc[idx[:cut]]; te=df.iloc[idx[cut:]]
    feature_sets={
        'D_only':['compression_deficit'],
        'D_plus_residual':['compression_deficit','max_residual'],
        'full_controls':['compression_deficit','max_residual','total_stopping_time','odd_steps','log_n']}
    rows=[]; full=None
    for name,features in feature_sets.items():
        beta,r2tr,r2te,means,stds,pred=fit_ols(tr,te,features)
        row={'model':name,'features':';'.join(features),'train_r2':r2tr,'test_r2':r2te,'intercept':beta[0]}
        for f,b in zip(features,beta[1:]): row['std_beta_'+f]=b
        rows.append(row)
        if name=='full_controls': full=(beta,r2tr,r2te,features,pred,te)
    mdf=pd.DataFrame(rows); mdf.to_csv(root/'tables'/'multivariate_holdout_models.csv',index=False)
    rel=pd.DataFrame([{'spearman_rho':rho,'spearman_p':p,'n':len(df),'holdout_fraction':.30,
                       'D_only_test_r2':float(mdf.loc[mdf.model=='D_only','test_r2'].iloc[0]),
                       'full_controls_test_r2':float(mdf.loc[mdf.model=='full_controls','test_r2'].iloc[0])}]); rel.to_csv(root/'tables'/'deficit_excursion_relation.csv',index=False)

    sample=df.sample(min(60000,len(df)),random_state=SEED)
    b=np.polyfit(df.compression_deficit,df.log_excursion_ratio,1)
    plt.figure(figsize=(7,4.8)); plt.scatter(sample.compression_deficit,sample.log_excursion_ratio,s=3,alpha=.18)
    xs=np.linspace(sample.compression_deficit.min(),sample.compression_deficit.max(),100); plt.plot(xs,b[1]+b[0]*xs,linewidth=1.5)
    plt.xlabel('Maximum compression deficit'); plt.ylabel(r'$\log(M(n)/n)$'); plt.tight_layout(); plt.savefig(root/'figures'/'deficit_vs_excursion.pdf'); plt.close()

    beta,r2tr,r2te,features,pred,te=full
    plt.figure(figsize=(6,5)); plt.scatter(te.log_excursion_ratio,pred,s=3,alpha=.2); lo=min(te.log_excursion_ratio.min(),pred.min()); hi=max(te.log_excursion_ratio.max(),pred.max()); plt.plot([lo,hi],[lo,hi],linewidth=1)
    plt.xlabel('Observed holdout log excursion'); plt.ylabel('Predicted holdout log excursion'); plt.tight_layout(); plt.savefig(root/'figures'/'holdout_prediction.pdf'); plt.close()

    t=np.arange(0,int(df.total_stopping_time.max())+1); surv=np.array([(df.total_stopping_time>tt).mean() for tt in t]); pd.DataFrame({'t':t,'survival':surv}).to_csv(root/'tables'/'stopping_survival.csv',index=False)
    plt.figure(figsize=(7,4.5)); plt.semilogy(t,np.maximum(surv,1/len(df))); plt.xlabel('Raw Collatz steps t'); plt.ylabel(r'$P(\tau>t)$'); plt.tight_layout(); plt.savefig(root/'figures'/'stopping_survival.pdf'); plt.close()
    return df,summary,rel,mdf


def validation_ranges(n_each:int,root:Path):
    ranges=[10**6,10**9,10**12]
    rows=[]
    for base in ranges:
        df=trajectory_metrics_range(base,n_each)
        rho,_=spearmanr(df.compression_deficit,df.log_excursion_ratio)
        rng=np.random.default_rng(SEED+base%100000); idx=rng.permutation(len(df)); cut=int(.7*len(df)); tr=df.iloc[idx[:cut]]; te=df.iloc[idx[cut:]]
        features=['compression_deficit','max_residual','total_stopping_time','odd_steps','log_n']
        beta,r2tr,r2te,means,stds,pred=fit_ols(tr,te,features)
        rows.append({'range_start':base,'range_end':base+n_each-1,'N':n_each,'spearman_rho_D_Y':rho,'full_controls_test_r2':r2te,
                     'std_beta_D':beta[1],'std_beta_Rmax':beta[2],'std_beta_stopping':beta[3],'std_beta_odd_steps':beta[4],'std_beta_log_n':beta[5]})
    out=pd.DataFrame(rows); out.to_csv(root/'tables'/'out_of_sample_ranges.csv',index=False)
    plt.figure(figsize=(7,4.5)); x=np.arange(len(out)); plt.plot(x,out.spearman_rho_D_Y,marker='o',label='Spearman rho(D,Y)'); plt.plot(x,out.full_controls_test_r2,marker='s',label='Holdout R2 (full controls)'); plt.xticks(x,[r'$10^6$',r'$10^9$',r'$10^{12}$']); plt.xlabel('Range start'); plt.ylabel('Statistic'); plt.legend(); plt.tight_layout(); plt.savefig(root/'figures'/'cross_range_stability.pdf'); plt.close()
    return out


def residual_extremal_table(root:Path,n0:int=10**12+1):
    rows=[]
    for m in [8,16,32,64]:
        for offset in [0,1,2]:
            K=math.ceil(m*math.log(3,2))+offset
            lo,hi=residual_kernel_extrema(m,K); rlo,rhi=residual_bounds_fixed_total(n0,m,K); C=K*LN2-m*LN3
            if C>0:
                n_all_nondescend=lo/(3.0*(math.exp(C)-1.0))
                n_all_descend=hi/(3.0*(math.exp(C)-1.0))
            else:
                n_all_nondescend=math.inf; n_all_descend=math.inf
            rows.append({'m':m,'K':K,'K_over_m':K/m,'compression_C':C,'kernel_min':lo,'kernel_max':hi,'kernel_ratio':hi/lo,
                         'R_min_at_n0':rlo,'R_max_at_n0':rhi,'uniform_descent_certificate':C>rhi,
                         'n_below_all_words_nondescend_bound':n_all_nondescend,'n_above_all_words_descend_bound':n_all_descend})
    df=pd.DataFrame(rows); df.to_csv(root/'tables'/'sharp_residual_extrema.csv',index=False)
    plt.figure(figsize=(7,4.5))
    for off,g in df.groupby(df.groupby('m').cumcount()):
        plt.plot(g.m,np.log10(g.kernel_ratio),marker='o',label=f'K-ceil(m log2 3)={off}')
    plt.xlabel('Word length m'); plt.ylabel(r'$\log_{10}(S_{max}/S_{min})$'); plt.legend(); plt.tight_layout(); plt.savefig(root/'figures'/'residual_ordering_amplification.pdf'); plt.close()
    return df




def near_extremal_bound_table(root:Path, epsilon:float=0.20):
    """Export the theorem-generated near-extremal count bounds used in the manuscript."""
    rows=[]
    for m,K in [(8,13),(16,26),(32,51),(64,102)]:
        rows.append(near_extremal_count_bound(m,K,epsilon))
    df=pd.DataFrame(rows)
    df.to_csv(root/'tables'/'near_extremal_count_bounds.csv',index=False)
    return df

def arithmetic_residue_table(root:Path, spectrum_base:int=10**30+1, nstarts:int=200000):
    """Generate exact residue classes for the sharp front/back extremal words."""
    rows=[]
    for m in [8,16,32,64]:
        K=math.ceil(m*math.log(3,2))
        back,front=extremal_valuation_words(m,K)
        for name,word in [('front',front),('back',back)]:
            residue,modulus=valuation_word_residue(word)
            rows.append({
                'm':m,'K':K,'extremizer':name,'word':','.join(map(str,word)),
                'residue':str(residue),'modulus':str(modulus),
                'log10_modulus':math.log10(modulus),'normalized_residue':residue/modulus,
                'density_among_odds':valuation_word_density_among_odds(word),
                'expected_count_in_spectrum_window':nstarts*valuation_word_density_among_odds(word),
                'exact_count_in_spectrum_window':count_word_in_odd_progression(word, spectrum_base if spectrum_base%2 else spectrum_base+1, nstarts)
            })
    df=pd.DataFrame(rows)
    df.to_csv(root/'tables'/'extremal_residue_classes.csv',index=False)
    return df



def accelerated_paradoxical_prefix_census(nmax:int, root:Path, max_m:int=100):
    """Census positive-compression accelerated prefixes with n_m >= n_0.

    Starts are odd n in [3,nmax].  This is a finite diagnostic, not a
    convergence claim.  The exact residue theorem is used to verify every
    reported witness and to compute its word-level threshold.
    """
    events=[]
    pair_counts={}
    pair_starts={}
    for n0 in range(3, nmax+1, 2):
        x=n0; word=[]; K=0
        for m in range(1, max_m+1):
            x,k=accelerated_step(x); word.append(int(k)); K += int(k)
            if (1 << K) > 3**m and x >= n0:
                w=tuple(word)
                residue,modulus=valuation_word_residue(w)
                threshold=valuation_word_nondescend_threshold(w)
                count=valuation_word_witness_count(w)
                assert n0 % modulus == residue % modulus
                assert threshold is not None and n0 <= threshold
                events.append({'n0':n0,'m':m,'K':K,'n_m':x,'residue':str(residue),
                               'modulus':str(modulus),'threshold':float(threshold),
                               'witness_count_word':count,'witness_ratio':valuation_word_witness_ratio(w),
                               'word':','.join(map(str,w))})
                key=(m,K); pair_counts[key]=pair_counts.get(key,0)+1; pair_starts.setdefault(key,set()).add(n0)
            if x==1:
                break
    edf=pd.DataFrame(events)
    edf.to_csv(root/'data'/'accelerated_paradoxical_prefixes.csv.gz',index=False,compression='gzip')
    rows=[]
    for (m,K),cnt in sorted(pair_counts.items()):
        ss=pair_starts[(m,K)]
        rows.append({'m':m,'K':K,'events':cnt,'distinct_starts':len(ss),'min_start':min(ss),'max_start':max(ss)})
    pdf=pd.DataFrame(rows)
    pdf.to_csv(root/'tables'/'accelerated_paradoxical_pairs.csv',index=False)
    return edf,pdf

def generalized_phase(nodd:int,root:Path):
    starts=odd_starts(1,nodd); rows=[]
    for a in [1,3,5,7,9,11]:
        ks=[accelerated_step(s,a,1)[1] for s in starts]
        mean_k=float(np.mean(ks)); rows.append({'a':a,'mean_k_fresh_ensemble':mean_k,'mean_lambda_baseline':mean_k*LN2-math.log(a),'critical_k_log2a':math.log(a,2)})
    df=pd.DataFrame(rows); df.to_csv(root/'tables'/'generalized_phase.csv',index=False)
    plt.figure(figsize=(7,4.5)); plt.plot(df.a,df.mean_k_fresh_ensemble,marker='o',label='Observed fresh-start mean valuation'); plt.plot(df.a,df.critical_k_log2a,marker='s',label=r'Critical $\log_2 a$'); plt.xlabel('Odd multiplier a in an+1'); plt.ylabel('Valuation per odd step'); plt.legend(); plt.tight_layout(); plt.savefig(root/'figures'/'generalized_phase.pdf'); plt.close(); return df



def residue_threshold_sparsity_analysis(root:Path):
    """Generate theorem-level residue-threshold bounds and small exact checks."""
    rows=[]
    for m in [16,32,64,128]:
        for rho in [1.60,1.70,1.85,2.00]:
            K=max(m, math.ceil(rho*m))
            if (1 << K) <= 3**m:
                K=math.floor(m*math.log(3,2))+1
            d=witness_word_count_upper_bound(m,K)
            d['rho_actual']=K/m
            rows.append(d)
    df=pd.DataFrame(rows)
    df.to_csv(root/'tables'/'residue_threshold_word_bounds.csv',index=False)

    erows=[]
    for rho in [1.60,1.65,1.70,1.80,2.00,2.50]:
        if rho > math.log(3,2):
            erows.append(strict_supercritical_witness_exponent(rho))
    edf=pd.DataFrame(erows)
    edf.to_csv(root/'tables'/'residue_threshold_asymptotic_exponents.csv',index=False)

    exact=[]
    for m,K in [(3,5),(4,7),(5,8),(6,10),(7,12),(8,13),(8,14),(9,15)]:
        if (1 << K) > 3**m and math.comb(K-1,m-1) <= 500000:
            exact.append(exact_witness_word_count_small(m,K,max_words=500000))
    xdf=pd.DataFrame(exact)
    xdf.to_csv(root/'tables'/'residue_threshold_exact_small.csv',index=False)

    plt.figure(figsize=(7,4.8))
    for rho,g in df.groupby('rho_actual'):
        pass
    for target in [1.60,1.70,1.85,2.00]:
        gg=[]
        for m in [16,32,64,128]:
            K=max(m,math.ceil(target*m))
            if (1<<K)<=3**m:
                K=math.floor(m*math.log(3,2))+1
            b=witness_word_count_upper_bound(m,K)
            gg.append((m,max(b['fraction_bound'],1e-300)))
        plt.plot([x[0] for x in gg],[x[1] for x in gg],marker='o',label=fr'$\rho\approx{target:.2f}$')
    plt.yscale('log'); plt.xlabel('Word length m'); plt.ylabel('Rigorous upper bound on witness-word fraction')
    plt.legend(); plt.tight_layout(); plt.savefig(root/'figures'/'residue_threshold_sparsity.pdf'); plt.close()
    return df,edf,xdf


def main():
    ap=argparse.ArgumentParser(description='Reproduce all analyses for the Collatz compression manuscript.')
    ap.add_argument('--out',default='../analysis_output')
    ap.add_argument('--n-raw',type=int,default=200000)
    ap.add_argument('--n-spectrum',type=int,default=200000)
    ap.add_argument('--max-m',type=int,default=128)
    ap.add_argument('--n-validation',type=int,default=10000)
    ap.add_argument('--workers',type=int,default=min(8, os.cpu_count() or 1))
    args=ap.parse_args()
    root=Path(args.out).resolve(); ensure_dirs(root)
    canonical=[1,2,4,8,16,32,64,96,128]; ms=[m for m in canonical if m<=args.max_m]
    spectrum_analysis(args.n_spectrum,ms,args.max_m,root,workers=args.workers)
    trajectory_analysis(args.n_raw,root)
    validation_ranges(args.n_validation,root)
    residual_extremal_table(root)
    near_extremal_bound_table(root)
    residue_threshold_sparsity_analysis(root)
    arithmetic_residue_table(root, nstarts=args.n_spectrum)
    accelerated_paradoxical_prefix_census(args.n_raw, root)
    generalized_phase(min(args.n_spectrum,100000),root)
    print(root)

if __name__=='__main__': main()
