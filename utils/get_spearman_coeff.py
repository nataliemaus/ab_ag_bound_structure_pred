from scipy import stats


def get_spearman_r(x,y):
    res = stats.spearmanr(x, y)
    return res.statistic