from scipy import stats


def get_spearman_r(x,y):
    res = stats.spearmanr(x, y)
    try:
        r_val = res.statistic.item()
    except:
        try:
            r_val = res.statistic # is item wrong? 
        except:
            r_val = res.correlation # older v of scipy 
    return r_val 