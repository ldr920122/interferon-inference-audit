import numpy as np, pandas as pd
def ridge(X,Y,seed,alpha=10000,nrand=1000):
    rng=np.random.default_rng(seed)
    X=(X-X.mean())/X.std(); Y=(Y-Y.mean())/Y.std(); Xa=X.values
    P=np.linalg.solve(Xa.T@Xa+alpha*np.eye(Xa.shape[1]),Xa.T); b=P@Y.values
    s1=np.zeros_like(b); s2=np.zeros_like(b); Yv=Y.values
    for _ in range(nrand):
        br=P@Yv[rng.permutation(len(Yv))]; s1+=br; s2+=br**2
    mu=s1/nrand; sd=np.sqrt(s2/nrand-mu**2)
    return pd.DataFrame((b-mu)/sd,index=X.columns,columns=Y.columns)
