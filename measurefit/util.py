from time import time
from math import pi, acos, atan, sin, cos
import numpy as np
from numpy.linalg import inv

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

# 度分秒转弧度
def dms2r(v):
    return (v//1 + v%1*100//1/60 + v*10000%100/60/60)/180*pi

# 将度分秒观测转换为弧度观测，以进行平差
def dms2rs(meas):
    rst = []
    for mea in meas:
        if mea[0]=='ang':
            rst.append(mea[:4] + (dms2r(mea[4]), mea[5]))
        elif mea[0]=='dir':
            rst.append(mea[:3] + (dms2r(mea[3]), mea[4]))
        elif mea[0]=='ts':
            rst.append(mea[:3] + (dms2r(mea[3]), dms2r(mea[4]), mea[5]))
        else: rst.append(mea)
    return rst

# 气象改正
def count_ppm(t, h, p):
    # p *= 10
    x = (7.5 * t) / (237.3 + t) + 0.07857
    a = 1/ 273.16
    d = (0.29525 * p) / (1 + a * t)
    e = (4.126 * 1e-4 * h) / (1 + a * t) * 10 ** x
    return - (286.34 - d + e)

# 全站仪加乘常数，ppm改正
def ts_adjust(mea, b=0, k=1, ppm=0):
    if mea[0] != 'ts': return mea
    return *mea[:5], mea[5]*(1+ppm/1e6)*k+b/1000

# 批量全站仪加乘常数，ppm改正
def ts_adjusts(meas, b=0, k=1, ppm=0):
    return [ts_adjust(i, b, k, ppm) for i in meas]

# 全站仪转水平角，高差，平距,K是遮光系数，0表示仅作曲率改正，1表示不做任何改正
def ts2ahd(mea, db, dk, da, K=0):
    if mea[0] != 'ts': return mea
    _, p1, p2, ax, ay, l = mea[:6]

    d = l*cos(ay) - (1-K)/4/6369000*l**2*sin(2*ay) # 平距
    h = l*sin(ay) + (1-K)*d**2 /2/6369000 # 高差
    
    dd = (db**2+(l/1000*dk)**2)*cos(ay)**2
    dd += (da/60/60/180*pi * l*1000 * sin(ay))**2

    dh = (db**2+(l/1000*dk)**2)*sin(ay)**2
    dh += (da/60/60/180*pi * l*1000 * cos(ay))**2
    return ('axhd', p1, p2, ax, h, d, da, dh**0.5, dd**0.5)
    
# 批量全站仪平距，高差，折光差改正
def ts2ahds(meas, db, dk, da, K=0):
    return [ts2ahd(i, db, dk, da, K) for i in meas]

# 平距进行统一高程投影
def distproj(mea, h, h0):
    k = 1-(h-h0)/6369000
    return *mea[:5], mea[5]*k, *mea[6:]

# 批量平距统一高程投影
def distprojs(meas, h0):
    rst = []; h = None
    for mea in meas:
        if mea[0]=='station': h=mea[2] # 站点高程
        if mea[0]=='axhd':
            rst.append(distproj(mea, h + mea[4]/2, h0))
        else: rst.append(mea)
    return rst

# 将全站仪观测分解成角度，距离，高差观测
def splitmeas(meas):
    rst = []
    for i in range(len(meas)):
        m = meas[i]
        if not m[0] in ('dirl','ts', 'axhd'):
            rst.append(m); continue
        if m[0]=='dirl':
            rst.append(('dir', m[1], m[2], m[3], m[5]))
            rst.append(('dist', m[1], m[2], m[4], m[6]))
        if m[0]=='axhd':
            rst.append(('dir', m[1], m[2], m[3], m[6]))
            rst.append(('dist', m[1], m[2], m[5], m[8]))
            rst.append(('lev',m[1], m[2], m[4], m[7]))
    return rst
            
# 距离投影到水准高程
def dist_adjust(mea, h, h0=None):
    if h0 is None: return mea
    d = mea[4] * cos(mea[3]) * (1-(h-h0)/6369000)
    return (*mea[:3], 0, d)

# 高差改正
def level_adjust(mea, h, h0=None):
    if h0 is None:
        return (*mea[:2], mea[4] * sin(mea[3]), *mea[3:])
    d = mea[4] * cos(mea[3]) * (1-(h-h0)/6369000)
    h = mea[4] * sin(mea[3]) + d**2 / 2 / 6369000
    return (*mea[:2], h, *mea[3:])

# 针对2x2矩阵的快速求逆
def inv2x2(m):
    (a, b), (c, d) = m[0], m[1]
    det = a * d - b * c
    return np.array([[d/det, -b/det], [-c/det, a/det]])

# 针对3x3矩阵的快速求逆
def inv3x3(m):
    (a, b, c), (d, e, f), (g, h, i) = m[0], m[1], m[2]
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    return np.array([
        [(e * i - f * h)/det, (c * h - b * i)/det, (b * f - c * e)/det],
        [(f * g - d * i)/det, (a * i - c * g)/det, (c * d - a * f)/det],
        [(d * h - e * g)/det, (b * g - a * h)/det, (a * e - b * d)/det]
    ])

# 构建后方交会严密平差矩阵
def build_backward(pts, meas, d1=1, d2=1e-3, d3=1):
    pts = dict([(i[0], i) for i in pts])
    num_meas = len(meas)
    B = np.zeros((2 * num_meas, 3))
    l = np.zeros(2 * num_meas)
    l0 = np.zeros(2 * num_meas)
    P = np.zeros(2 * num_meas)

    for i, m in enumerate(meas):
        p1, p0 = pts[m[2]][1:4], pts[m[1]][1:4]
        a0 = p0[2]
        dv = minus(p1, p0)
        nv, nv2 = norm(dv), norm(dv[:2])
        kx, ky = dv[0] / nv2, dv[1] / nv2
        B[2 * i, :2] = [-kx, -ky]
        l[2 * i] = m[5] * 1000
        l0[2 * i] = nv2 * 1000
        P[2 * i] = d1**2 + (nv2 * d2)**2

        k = 180 * 60 * 60 / np.pi
        kx = -dv[1] / nv2**2 * k / 1000
        ky = dv[0] / nv2**2 * k / 1000
        B[2*i+1] = [-kx, -ky, -1]
        l[2*i+1] = ((m[3] + a0) % (pi*2)) * k
        r = angleX(dv)
        l0[2*i+1] = fita(r, (m[3] + a0) % (pi*2)) * k
        P[2*i+1] = d3 ** 2
    return B, l, l0, P


# 后方交汇求解
def backward_x(pts, meas, db=1, dk=1e-3, da=1, accu=1e-1, maxiter=3):
    meas = meas[1:] # 去掉station
    pts = [list(i) for i in pts]
    dic = dict([(i[0], i) for i in pts])
    p1, p0 = dic[meas[0][2]][1:4], dic[meas[0][1]][1:4]
    dic[meas[0][1]][3] = fita(angleX(minus(p1, p0)) - meas[0][3], 0)
        
    for i in range(maxiter):
        B, l, l0, p = build_backward(pts, meas, db, dk, da)
        # dx, Dx = count(B, l-l0, np.diag(1.0/p))
        Dx = inv3x3(B.T / p @ B)
        dx = Dx @ B.T / p @ (l-l0)
        pts[-1][1] += dx[0]/1000
        pts[-1][2] += dx[1]/1000
        pts[-1][3] += dx[2] * np.pi/180/60/60
        
        if np.abs(dx).max()<accu:
            L = B @ dx + l0 # 观测平差值
            sta = (l - L).T * 1/p @ (l - L) # 后验总误差
            h = [dic[m[2]][3]-m[4] for m in meas]
            dh = [1/m[7]**2 for m in meas]
            h0 = sum([i*j for i,j in zip(h,dh)])/sum(dh)
            
            sta += sum([(i-h0)**2*j for i,j in zip(h, dh)])*1e6
            if len(meas)==1: sta = 1
            else: sta = sta/(len(meas)*3-4)
            return (*pts[-1][:3], h0, False, *Dx[:2,:2].ravel(),
                    1/sum(dh)), sta**0.5
    return None

# 单点极坐标求解
def polar(pts, meas):
    pts = dict([(i[0], i) for i in pts])
    p0, p1 = meas[1][1], meas[1][2]
    p0, p1 = pts[p0], pts[p1]
    
    a0 = angleX(minus(p1[1:4], p0[1:4])) - meas[1][3]
    rst = []
    for mea in meas[2:]:
        ax, h, d = mea[3]+a0, mea[4], mea[5]
        x = p0[1] + d * cos(ax)
        y = p0[2] + d * sin(ax)
        z = p0[3] + h
        da, dh, dl = mea[6], mea[7]**2, mea[8]**2
        da = (da/60/60/180*pi*d*1000)**2*2
        
        m = np.array([[cos(ax), -sin(ax)],[sin(ax), cos(ax)]])
        Dx = m @ np.array([[dl, 0],[0, da]]) @ m.T
        rst.append((mea[2], x, y, z, False, *Dx.ravel(), dh))
    return rst

# 将多组极坐标结果融合成交汇成果
def merge_xyh(xyhs):
    if len(xyhs)==0: return xyhs[0]
    x0 = np.array([xyh[1:4] for xyh in xyhs])
    xs = (x0 - x0.mean(axis=0))*1000
    ps = [[[xyh[5], xyh[6], 0],
           [xyh[7], xyh[8], 0],
           [0,   0,   xyh[9]]] for xyh in xyhs]
    ps = np.array([inv3x3(i) for i in ps])
    dx = inv3x3(sum(ps))
    x = dx @ (ps @ xs[:,:,None]).sum(axis=0)
    x0 = x0.mean(axis=0) + x.ravel()/1000
    e = xs - x.ravel()
    sta = (e[:,None,:] @ ps @ e[:,:,None]).sum()
    if len(xyhs)==1: stafter = 1
    else: stafter = sta/(len(xyhs)*3-3)
    return (xyhs[0][0], *x0, False, 
            *dx[:2,:2].ravel(), dx[2,2]), stafter**0.5

# 批量极坐标融合成交汇成果
def merge_xyhs(xyhss):
    group = {}
    for xyh in [j for i in xyhss for j in i]:
        if not xyh[0] in group: group[xyh[0]] = []
        group[xyh[0]].append(xyh)
    return [merge_xyh(i)[0] for i in group.values()]

# 前方交会，观测列表需要符合多个极坐标拼接起来的形式
def forward_x(pts, meas):
    sep = [i for i in range(len(meas)) if meas[i][0]=='station']+[None]
    xyhs = [polar(pts, meas[s:e]) for s,e in zip(sep[:-1], sep[1:])]
    return merge_xyhs(xyhs)


if __name__ == '__main__':
    # 后方交汇
    points = [
        ('TN2', 1000, 1000, 101, True),
        ('TN8', 1367.77286, 1117.43005, 102, True),
        ('TN9', 1294.76794, 1789.72043, 103, True),
        ('TN1', 893, 1603, 100, False)]

    measures = [
        ('station', 'TN1', None),
        ('ts', 'TN1', 'TN2', 0, 0.00163, 612.41862),
        ('ts', 'TN1', 'TN8', 0.598815, 0.00294, 678.73168),
        ('ts', 'TN1', 'TN9', 1.8317007, 0.00678, 442.39598)
    ]
    
    measx = ts2ahds(measures, 1, 1, 1, K=1) # 全站仪数据转水平角，高差，平距
    rst, sta = backward_x(points, measx)
    
    points = [
        ('P1', 0, 0, 0, True),
        ('P2', 1000, 0, 0, True),
        ('Px', 1000, 577, 100, False)
    ]

    # A站观测
    meas_a = [
        ('station', 'P1', None),
        ('ts', 'P1', 'P2', 0, 0, 1000),
        ('ts', 'P1', 'Px', 0.52359, 0.0864, 1158.84, None),
        # 这里可以有更多未知点观测记录
    ]
    # B站观测
    meas_b = [
        ('station', 'P2', None),
        ('ts', 'P2', 'P1', 3.141592, 0, 1000, None),
        ('ts', 'P2', 'Px', 1.57079, 0.17151, 585.94, None),
        # 这里可以有更多未知点观测记录
    ]

    # 改正，转高差，平距
    meas = ts_adjusts(meas_a, b=0, k=1, ppm=0)
    meas = ts2ahds(meas, 1, 1, 1, K=1)
    xyh_a = polar(points, meas)

    meas = ts_adjusts(meas_b, b=0, k=1, ppm=0)
    meas = ts2ahds(meas, 1, 1, 1, K=1)
    xyh_b = polar(points, meas)

    xyh_ab = merge_xyhs([xyh_a, xyh_b])


    meas = ts_adjusts(meas_a+meas_b, b=0, k=1, ppm=0) # 加乘常数改正，气象改正
    meas = ts2ahds(meas, 1, 1, 1, K=1) # 全站仪数据转水平角，高差，平距
    
    xyh_abx = forward_x(points, meas)

