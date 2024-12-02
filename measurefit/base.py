import numpy as np
inv = np.linalg.inv
from math import pi, sin, cos, asin, acos

# 矢量减法
def minus(v1, v2):
    return v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]

# 矢量叉乘
def cross(v1, v2):
    return v1[0] * v2[1] - v2[0] * v1[1]

# 矢量点乘
def dot(v1, v2):
    return v1[0]*v2[0] + v1[1]*v2[1]

# 计算方位角
def angleX(v):
    a = acos(v[0] / (v[0]*v[0] + v[1]*v[1])**0.5)
    return a if v[1]>0 else pi * 2 - a

# 矢量求模
def norm(v):
    return (v[0]*v[0] + v[1]*v[1])**0.5

# 角度临界值匹配
def fita(a1, a2):
    if a1 - a2 > pi: a1 -= pi * 2
    if a2 - a1 > pi: a1 += pi * 2
    return a1

# 距离观测
def distance(p1, p2, dis, sta, db, dk, lw):
    dv = minus(p2[1:4], p1[1:4])
    l, limit = norm(dv), {}
    if not p2[4]:
        limit[p2[0]+'.x'] = dv[0]/l
        limit[p2[0]+'.y'] = dv[1]/l
    if not p1[4]:
        limit[p1[0]+'.x'] = -dv[0]/l
        limit[p1[0]+'.y'] = -dv[1]/l
    if sta is None: # 两种距离定权方式
        if lw==2: sta = db**2 + (dis/1000*dk)**2
        if lw==1: sta = (db + (dis/1000*dk))**2
    else: sta **= 2 # 给定精度
    return limit, l*1000, dis*1000, sta

# 方向观测
def direction(p1, p2, st, ang, sta, da):
    dv = minus(p2[1:4], p1[1:4])
    if st and st[1] is None: # 计算方位角
        st[1] = fita(angleX(dv)-ang, 0)
    l, limit = norm(dv), {}
    k = 180*60*60/1000/pi
    if not p2[4]:
        limit[p2[0]+'.x'] = -dv[1]/l**2*k
        limit[p2[0]+'.y'] = dv[0]/l**2*k
    if not p1[4]:
        limit[p1[0]+'.x'] = dv[1]/l**2*k
        limit[p1[0]+'.y'] = -dv[0]/l**2*k
    if sta!=0: limit[p1[0]+'.a0'] = -1
    k = 1 / pi*180*60*60
    sta = da**2 if sta is None else sta**2
    ang = (ang + (st[1] if st else 0))%(pi*2)
    ang0 = fita(angleX(dv), ang)
    return limit, ang0*k, ang*k, sta

# 角度观测
def angle(p1, p0, p2, ang, sta, da):
    dv1 = minus(p1[1:4], p0[1:4])
    dv2 = minus(p2[1:4], p0[1:4])
    l1, l2 = norm(dv1), norm(dv2)
    sign = (-1,1)[cross(dv1,dv2)>0]

    limit = {}
    k = 180*60*60/1000/pi
    if not p1[4]:
        limit[p1[0]+'.x'] = dv1[1]/l1**2*k*sign
        limit[p1[0]+'.y'] = -dv1[0]/l1**2*k*sign
    if not p0[4]:
        limit[p0[0]+'.x'] = -dv1[1]/l1**2*k*sign
        limit[p0[0]+'.x'] += dv2[1]/l2**2*k*sign
        limit[p0[0]+'.y'] = dv1[0]/l1**2*k*sign
        limit[p0[0]+'.y'] -= dv2[0]/l2**2*k*sign
    if not p2[4]:
        limit[p2[0]+'.x'] = -dv2[1]/l2**2*k*sign
        limit[p2[0]+'.y'] = dv2[0]/l2**2*k*sign
        
    k = 1 / pi*180*60*60
    sta = da**2 if sta is None else sta**2
    ang0 = acos(dot(dv1,dv2)/norm(dv1)/norm(dv2))
    return limit, ang0*k, ang*k, sta

# 水准观测
def level(p1, p2, dh, sta):
    limit = {}
    if not p1[4]: limit[p1[0]+'.h'] = -1
    if not p2[4]: limit[p2[0]+'.h'] = 1
    sta = -sta if sta<0 else sta**2
    return limit, (p2[3]-p1[3])*1000, dh*1000, sta

# 观测转约束表达式
def trans(pts, meas, db, dk, da, lw):
    t_meas, tags = [], []
    for mea in meas:
        if mea[0]=='dist' and not (pts[mea[1]][4] and pts[mea[2]][4]):
            para = pts[mea[1]], pts[mea[2]], mea[3], mea[4], db, dk, lw
            tags.append(['distance', mea[1], mea[2], mea[3]])
            t_meas.append(distance(*para))
        if mea[0]=='station' and not ('st:'+mea[1]) in pts:
            pts['st:'+mea[1]] = ['st:'+mea[1], None]
        if mea[0]=='dir':
            st = 'st:'+mea[1] if mea[4]!=0 else None
            para = pts[mea[1]], pts[mea[2]], pts.get(st), mea[3], mea[4], da
            tags.append(['direction', mea[1], mea[2], mea[3]/pi*180])
            t_meas.append(direction(*para))
        if mea[0]=='lev' and not (pts[mea[1]][4] and pts[mea[2]][4]):
            tags.append(['level', mea[1], mea[2], mea[3]])
            t_meas.append(level(pts[mea[1]], pts[mea[2]], mea[3], mea[4]))
        if mea[0]=='ang':
            tags.append(['angle', mea[1]+'-'+mea[2]+'-'+mea[3], mea[4]])
            para = pts[mea[1]], pts[mea[2]], pts[mea[3]], mea[4], mea[5], da
            t_meas.append(angle(*para))
    return t_meas, tags

# 更新概略坐标
def update(pts, dx, idx):
    for p_i, v in zip(idx, dx):
        p, i = p_i.split('.')
        if i=='a0': pts['st:'+p][1] += v*pi/(180*60*60)
        else: pts[p]['nxyh'.index(i)] += v/1000
        
# 构建矩阵
def build(meas, idx):
    B = np.zeros((len(meas), len(idx)))
    l, l0, p = np.zeros((3, len(meas)))
    for n, mea in enumerate(meas):
        for k,v in mea[0].items():
            B[n, idx.index(k)] = v
        l0[n], l[n], p[n] = mea[1], mea[2], mea[3]
    return B, l0, l, p

# Bx=l:p最小二乘解
def count(B, l, p):
    Dx = inv(B.T @ p @ B)
    x = Dx @ B.T @ p @ l
    return x, Dx

# Cx=b 约束下的 Bx=l:p最小二乘
def count_lim(B, l, p, C, b):
    Dx = inv(B.T @ p @ B)
    Dc = inv(C @ Dx @ C.T)
    x = Dx @ B.T @ p @ l
    r = Dc @ (C @ x - b)
    x = x - Dx @ C.T @ r
    Dx -= Dx @ C.T @ Dc @ C @ Dx.T
    return x, Dx
    
# 单步求解
def accept(x, dx, e, m, b, d):
    dx, d = inv(dx), inv(d)
    mtd = m.T @ d
    newd = inv(dx + mtd @ m)
    newx = newd @ (dx @ x + mtd @ b)
    # e1 = newx - x
    # e1 = e1.T @ dx @ e1
    # e2 = m @ newx - b
    # e2 = e2.T @ d @ e2
    return newx, newd #, e+e1+e2

# 递推求解
def icount(B, l, p):
    x, e = np.zeros(B.shape[1]), 0
    dx = np.diag(np.ones(len(x)))*1e6
    para = B[:,None], l[:,None], p[:,None,None]
    for m, b, d in zip(*para):
        x, dx = accept(x, dx, e, m, b, d)
    return x, dx

# 求解函数
def solve(pts, meas, maxiter=3, accu=0.1, db=1, dk=1, da=1, iter=False, lw=2):
    pts = dict((i[0], list(i)) for i in pts)
    for i in range(maxiter):
        t_meas, tags = trans(pts, meas, db, dk, da, lw)
        idx = [list(i[0].keys()) for i in t_meas]
        idx = sorted(set(sum(idx, [])))
        
        B, l0, l, p = build(t_meas, idx)
        msk = p==0;
        if msk.sum()==0:
            if iter: x, Dx = icount(B, l-l0, p)
            else: x, Dx = count(B, l-l0, np.diag(1/p))
        else: # 含有约束条件
            C, b0, b = B[msk], l0[msk], l[msk]
            B, l0, l, p = B[~msk], l0[~msk], l[~msk], p[~msk]
            x, Dx = count_lim(B, l-l0, np.diag(1/p), C, b-b0)
        # elif iter: x, Dx = icount(B, l-l0, p)
        update(pts, x, idx)

        # 求解成功，组装成果
        if np.abs(x).max()<accu:
            L = B @ x + l0 # 观测平差值
            Dl = (B[:,None]@ Dx @B[:,:,None]).ravel() # 观测平差精度
            stafter = (l - L).T * 1/p @ (l - L) # 后验总误差
            if len(msk)-len(idx)==0: stafter = 1
            else: stafter /= len(msk)-len(idx) # 单位方差
            key = {'angle':3600, 'direction':3600, 'distance':1000, 'level':1000}
            for i, l, dl in zip(np.where(~msk)[0], L, Dl):
                tags[i].extend([l / key[tags[i][0]], dl])
            for i, l in zip(np.where(msk)[0], C@x+b0 if msk.sum() else []):
                tags[i].extend([l / key[tags[i][0]], 0])
            for p in pts.values():
                if p[0]+'.x' in idx:
                    ix, iy = [idx.index(p[0]+i) for i in ('.x', '.y')]
                    p.extend(Dx[[ix, ix, iy, iy], [ix, iy, ix, iy]])
                else: p.extend([0, 0, 0, 0])
                if p[0]+'.h' in idx:
                    ih = idx.index(p[0]+'.h')
                    p.append(Dx[ih, ih])
                else: p.append(0)
                
            pts = [tuple(i) for i in pts.values() if not 'st:' in i[0]]
            return pts, sorted([tuple(i) for i in tags]), stafter
    raise Exception('solve failed')

if __name__ == '__main__':
    pts =[('TN2', 1000, 1000, 0, True),
          ('TN8', 1367.77286, 1117.43005, 0, True),
          ('TN9', 1294.76794, 1789.72043, 0, True),
          ('TN1', 893, 1603, 0, False)]

    meas = [('station', 'TN1', None),
          ('dir', 'TN1', 'TN2', 0, None),
          ('dist', 'TN1', 'TN2', 612.41862, None),
          ('dir', 'TN1', 'TN8', (34+18/60+34.46/60/60)/180*pi, None),
          ('dist', 'TN1', 'TN8', 678.73168, None),
          ('dir', 'TN1', 'TN9', (104+56/60+55.39/60/60)/180*pi, None),
          ('dist', 'TN1', 'TN9', 442.39598, None)]

    from time import time
    start = time()
    for i in range(100):
        rst, info, sta = solve(pts, meas, maxiter=3, accu=0.1, db=1, dk=1, da=1, iter=False)
    print(time()-start)
