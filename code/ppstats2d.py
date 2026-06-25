"""2D point-pattern diagnostics for the hyperuniformity-class characterization.
Data-agnostic; validated on synthetic patterns before use on real fields."""
import numpy as np
from scipy.spatial import cKDTree

# ---------- nulls / synthetic generators ----------
def csr(n, L, rng):              # complete spatial randomness (Poisson)
    return rng.uniform(0, L, size=(n, 2))

def rsa_hardcore(n, L, dmin, rng, maxtry=200000):
    """random sequential adsorption hard disks (liquid-like, repulsive)."""
    pts=[]; tree=None; t=0
    while len(pts)<n and t<maxtry:
        p=rng.uniform(0,L,size=2); t+=1
        if not pts or np.min(np.linalg.norm(np.array(pts)-p,axis=1))>=dmin:
            pts.append(p)
    return np.array(pts)

def perturbed_lattice(m, L, sigma, rng):
    """m x m square lattice + Gaussian jitter -> disordered hyperuniform."""
    a=L/m; g=(np.arange(m)+0.5)*a
    X,Y=np.meshgrid(g,g); P=np.column_stack([X.ravel(),Y.ravel()])
    return (P+rng.normal(0,sigma,P.shape))%L

# ---------- diagnostics ----------
def clark_evans(points, L):
    """NNI R = mean NN dist / expected(CSR). R<1 cluster, =1 random, >1 regular."""
    n=len(points); tree=cKDTree(points)
    d,_=tree.query(points,k=2); dnn=d[:,1]
    rho=n/(L*L); expected=1.0/(2*np.sqrt(rho))
    R=dnn.mean()/expected
    z=(dnn.mean()-expected)/(0.26136/np.sqrt(n*rho))
    return R, z, dnn.mean()

def pair_correlation(points, L, rmax, dr):
    """g(r) via radial histogram with periodic-free edge handling (density norm)."""
    n=len(points); rho=n/(L*L); tree=cKDTree(points)
    edges=np.arange(dr,rmax+dr,dr); g=np.zeros(len(edges)-1)
    centers=0.5*(edges[1:]+edges[:-1])
    counts=np.zeros(len(edges)-1)
    for i,p in enumerate(points):
        idx=tree.query_ball_point(p,rmax)
        d=np.linalg.norm(points[idx]-p,axis=1); d=d[d>0]
        h,_=np.histogram(d,bins=edges); counts+=h
    ring=np.pi*(edges[1:]**2-edges[:-1]**2)
    g=counts/(n*rho*ring)
    return centers, g

def structure_factor(points, L, kmax, nk=60):
    """radially-averaged S(k)=<|sum exp(-i k.r)|^2>/N on allowed k=2pi n/L grid."""
    n=len(points)
    nmax=int(np.ceil(kmax*L/(2*np.pi)))
    ns=np.arange(-nmax,nmax+1)
    KX,KY=np.meshgrid(ns,ns); 
    kx=2*np.pi*KX.ravel()/L; ky=2*np.pi*KY.ravel()/L
    kmag=np.sqrt(kx**2+ky**2)
    sel=(kmag>0)&(kmag<=kmax)
    kx,ky,kmag=kx[sel],ky[sel],kmag[sel]
    # S(k) per allowed vector
    ph=points[:,0][None,:]*kx[:,None]+points[:,1][None,:]*ky[:,None]
    S=(np.abs(np.exp(-1j*ph).sum(axis=1))**2)/n
    # radial average
    kb=np.linspace(0,kmax,nk+1); kc=0.5*(kb[1:]+kb[:-1]); Sr=np.full(nk,np.nan)
    for i in range(nk):
        m=(kmag>=kb[i])&(kmag<kb[i+1])
        if m.sum()>0: Sr[i]=S[m].mean()
    return kc, Sr

def number_variance(points, L, radii, nsamp, rng):
    """sigma^2(R) = var of count in random disks radius R (interior, no edge bias)."""
    tree=cKDTree(points); out=[]
    for R in radii:
        c=rng.uniform(R,L-R,size=(nsamp,2))
        cnt=np.array([len(tree.query_ball_point(p,R)) for p in c])
        out.append(cnt.var())
    return np.array(radii), np.array(out)

# ---------- rectangular-window versions for real fields ----------
def structure_factor_rect(points, Lx, Ly, kmax, nk=40):
    n=len(points)
    nxmax=int(np.ceil(kmax*Lx/(2*np.pi))); nymax=int(np.ceil(kmax*Ly/(2*np.pi)))
    nx=np.arange(-nxmax,nxmax+1); ny=np.arange(-nymax,nymax+1)
    KX,KY=np.meshgrid(nx,ny); kx=2*np.pi*KX.ravel()/Lx; ky=2*np.pi*KY.ravel()/Ly
    km=np.sqrt(kx**2+ky**2); sel=(km>0)&(km<=kmax); kx,ky,km=kx[sel],ky[sel],km[sel]
    ph=points[:,0][None,:]*kx[:,None]+points[:,1][None,:]*ky[:,None]
    S=(np.abs(np.exp(-1j*ph).sum(axis=1))**2)/n
    kb=np.linspace(0,kmax,nk+1); kc=0.5*(kb[1:]+kb[:-1]); Sr=np.full(nk,np.nan)
    for i in range(nk):
        m=(km>=kb[i])&(km<kb[i+1])
        if m.sum()>0: Sr[i]=S[m].mean()
    return kc,Sr

def number_variance_rect(points, box, radii, nsamp, rng):
    xmin,xmax,ymin,ymax=box; tree=cKDTree(points); out=[]
    for R in radii:
        if xmax-xmin<2*R or ymax-ymin<2*R: out.append(np.nan); continue
        c=np.column_stack([rng.uniform(xmin+R,xmax-R,nsamp),rng.uniform(ymin+R,ymax-R,nsamp)])
        cnt=np.array([len(tree.query_ball_point(p,R)) for p in c]); out.append(cnt.var())
    return np.array(radii),np.array(out)

def inhomogeneous_poisson(points, box, ncell, rng):
    """resample CSR within a coarse density grid matched to the data (heterogeneity null)."""
    xmin,xmax,ymin,ymax=box
    hx=np.linspace(xmin,xmax,ncell+1); hy=np.linspace(ymin,ymax,ncell+1)
    H,_,_=np.histogram2d(points[:,0],points[:,1],bins=[hx,hy])
    out=[]
    for ix in range(ncell):
        for iy in range(ncell):
            k=int(H[ix,iy])
            if k>0:
                out.append(np.column_stack([rng.uniform(hx[ix],hx[ix+1],k),rng.uniform(hy[iy],hy[iy+1],k)]))
    return np.vstack(out)
