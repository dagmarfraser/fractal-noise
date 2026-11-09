import sys, numpy as np, yasa, inspect, warnings
warnings.filterwarnings("ignore")
X = np.loadtxt(sys.argv[1], delimiter=",")[:, :20]
print(inspect.signature(yasa.irasa))
def run(**kw):
    out = []
    for r in range(X.shape[1]):
        _, _, _, fit = yasa.irasa(X[:, r], sf=100, band=(1, 20), return_fit=True, verbose=False, **kw)
        out.append(-fit["Slope"].iloc[0])
    return f"median {np.median(out):.3f} sd {np.std(out):.3f}"
print("default win4       ", run())
print("win10.24           ", run(win_sec=10.24))
print("win40.96           ", run(win_sec=40.96))
print("win10.24 detrend lin", run(win_sec=10.24, kwargs_welch=dict(average="median", window="hann", detrend="linear")))
print("win4 detrend lin   ", run(kwargs_welch=dict(average="median", window="hann", detrend="linear")))
print("win10.24 mean avg  ", run(win_sec=10.24, kwargs_welch=dict(average="mean", window="hann")))
