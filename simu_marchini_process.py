import numpy as np


def Xmat(nind, nsnp):
    return np.random.randint(0, 3, size=(nind, nsnp)).astype(np.float64)


def half_khatri_rao(x):
    out = np.zeros((x.shape[0], x.shape[1] * (x.shape[1] - 1) // 2))
    cnt = 0
    for i in range(x.shape[1] - 1):
        for j in range(i + 1, x.shape[1]):
            out[:, cnt] = x[:, i] * x[:, j]
            cnt += 1
    return out


def Intermat(X):
    inter = np.zeros((X.shape[0], X.shape[1] * (X.shape[1] + 1) // 2))
    inter[:, :nsnp] = X
    inter[:, nsnp:] = half_khatri_rao(X)

    return inter


def multiplicative_within_between(alpha, theta1, theta2, X1, X2):
    return alpha * ((1 + theta1) ** X1) * ((1 + theta2) ** X2)


def multiplicative_two_locus(alpha, theta, gamma, X1, X2):
    return alpha * (1 + theta) ** (X1 * X2)


def threshold_two_locus(alpha, theta, gamma, X1, X2):
    return alpha * (1 + theta) ** (
        np.asarray([1 if val > 1 else 0 for val in np.minimum(X1, X2)])
    )


def position(ii, jj, n):
    i = min(ii, jj)
    j = max(ii, jj)
    return i * n - i * (i + 1) // 2 + j - i - 1


def simulate(X, s, sp0, sp1, functions):
    nind, nsnp = X.shape
    Y = np.random.normal(0, s, (nind,))
    support_1d = np.arange(nsnp)
    np.random.shuffle(support_1d)
    support_1d = support_1d[:sp0]
    support_2d = []
    for i in range(sp1):
        indexes = np.arange(nsnp)
        np.random.shuffle(indexes)
        tmp = (min(indexes[0], indexes[1]), max(indexes[0], indexes[1]))
        if tmp in support_2d:
            i -= 1
        else:
            support_2d.append(tmp)
    if sp0 > 0:
        for idx in support_1d:
            val = (np.random.binomial(1, 0.5) * 2 - 1) * 0.1
            Y += X[:, idx] * val
    itt = []
    for idx in support_2d:
        val = (np.random.binomial(1, 0.5) * 2 - 1) * 0.1
        v1 = np.random.uniform(0.2, 0.5)
        v2 = np.random.uniform(0.2, 0.5)
        fun_idx = np.random.randint(0, len(functions))
        Y += functions[fun_idx](val, v1, v2, X[:, idx[0]], X[:, idx[1]])
        itt.append(fun_idx)

    Y -= np.mean(Y)
    Y /= np.std(Y)

    X = X / 2
    inter = Intermat(X)
    for col in range(inter.shape[1]):
        inter[:, col] -= np.mean(inter[:, col])

    sup2 = [position(v[0], v[1], nsnp) for v in support_2d]
    sup2f = []
    for i_f in range(len(functions)):
        sup2f.append([p for i, p in enumerate(sup2) if itt[i] == i_f])

    return inter, Y, support_1d, support_2d, sup2, sup2f


def save(X, inter, Y, sup1d, sup2d, sup2_by_fun, nsimu, nind, spv, sp, s, funs):
    # create directory with parameters as name
    import os
    import json

    dir = "results"
    if not os.path.exists(dir):
        os.mkdir(dir)
    dir = dir + "/" + f"{nind}_{sp[0]}_{sp[1]}_{spv}_{s}_{funs}_{nsimu}"
    if not os.path.exists(dir):
        os.mkdir(dir)

    np.save(dir + "/X", X)
    np.save(dir + "/Y", Y)
    np.save(dir + "/inter", inter)
    mydict = {
        "sup1d": [int(i) for i in sup1d],
        "sup2d": [(int(i), int(j)) for i, j in sup2d],
    }
    if funs == "0":
        mydict.update({"sup2d_t1": [int(k) for k in sup2_by_fun[0]]})
    elif funs == "1":
        mydict.update({"sup2d_t2": [int(k) for k in sup2_by_fun[0]]})
    elif funs == "2":
        mydict.update({"sup2d_t3": [int(k) for k in sup2_by_fun[0]]})
    else:
        mydict.update(
            {
                "sup2d_t1": [int(k) for k in sup2_by_fun[0]],
                "sup2d_t2": [int(k) for k in sup2_by_fun[1]],
                "sup2d_t3": [int(k) for k in sup2_by_fun[2]],
            }
        )

    with open(dir + "/meta.json", "w") as f:
        json.dump(mydict, f)


if __name__ == "__main__":
    nsnp = 10
    nsimu = 3
    functions = [
        multiplicative_within_between,
        multiplicative_two_locus,
        threshold_two_locus,
    ]

    NIND = {"low accession count": 10, "high accession count": 20}
    SPV = {"low signal count": 1, "high signal count": 10}
    SP = {"with 1D signal": (2, 3), "without 1D signal": (0, 5)}
    S = {"low noise": 0.1, "high noise": 5}

    for nind in NIND.items():
        for spv in SPV.items():
            for sp in SP.items():
                for s in S.items():
                    for i in range(nsimu):
                        for ii in range(len(functions)):
                            X = Xmat(nind[1], nsnp)
                            inter, Y, support_1d, support_2d, sup2, sup2f = simulate(
                                X,
                                s[1],
                                sp[1][0] * spv[1],
                                sp[1][1] * spv[1],
                                [functions[ii]],
                            )
                            save(
                                X,
                                inter,
                                Y,
                                support_1d,
                                support_2d,
                                sup2f,
                                i,
                                nind[1],
                                spv[1],
                                sp[1],
                                s[1],
                                f"{ii}",
                            )

                        X = Xmat(nind[1], nsnp)
                        inter, Y, support_1d, support_2d, sup2, sup2f = simulate(
                            X, s[1], sp[1][0] * spv[1], sp[1][1] * spv[1], functions
                        )
                        save(
                            X,
                            inter,
                            Y,
                            support_1d,
                            support_2d,
                            sup2f,
                            i,
                            nind[1],
                            spv[1],
                            sp[1],
                            s[1],
                            "-1",
                        )
