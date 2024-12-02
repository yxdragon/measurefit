# 转 dict
def mapping(pts=None, info=None, sta=0.0):
    if not pts is None:
        key = ['name', 'x', 'y', 'h', 'type', 'dxx', 'dxy', 'dyx', 'dyy', 'dhh']
        pts = [dict(zip(key, i)) for i in pts]
    if not info is None:
        key1 = ['measure', 'from', 'to', 'l', 'dl']
        key2 = ['measure', 'from-to', 'l', 'dl']
        info = [dict(zip(key2 if i[0]=='angle' else key1, i)) for i in info]
    if pts is None and not info is None: return info
    if info is None and not pts is None: return pts
    return pts, info
    
# 生成报告
def report(pts, info=None, sta=0.0):
    cont = []
    cont.append('点位信息及精度:')
    key = ['name', 'x', 'y', 'h', 'type', 'dxx', 'dxy', 'dyy', 'dhh']
    temp = '%-8s %-12s %-12s %-12s %-8s %-10s %-10s %-10s %-10s'
    cont.append(temp%tuple(key)); cont.append('-'*98)
    temp = '%-8s %-12.4f %-12.4f %-12.4f %-8s %-10.4f %-10.4f %-10.4f %-10.4f'
    for item in pts: cont.append(temp%tuple(item[:-3]+item[-2:]))

    if info and 'angle' in [i[0] for i in info]:
        cont.append('\n角度观测及精度:')
        key = ['measure', 'from_to', 'mea(°)', 'value(°)', 'dl(″)']
        temp = '%-10s %-8s %-12s %-12s %-12s'
        cont.append(temp%tuple(key)); cont.append('-'*63)
        temp = '%-10s %-8s %-12.4f %-12.4f %-12.4f'
        for item in [i for i in info if i[0]=='angle']:
            cont.append(temp%tuple(item))
            
    if info and 'direction' in [i[0] for i in info]:
        cont.append('\n方向观测及精度:')
        key = ['measure', 'from', 'to', 'mea(°)', 'value(°)', 'dl(″)']
        temp = '%-10s %-8s %-8s %-12s %-12s %-12s'
        cont.append(temp%tuple(key)); cont.append('-'*63)
        temp = '%-10s %-8s %-8s %-12.4f %-12.4f %-12.4f'
        for item in [i for i in info if i[0]=='direction']:
            cont.append(temp%tuple(item))

    if info and 'distance' in [i[0] for i in info]:
        cont.append('\n边长观测及精度:')
        key = ['measure', 'from', 'to', 'mea(m)', 'value(m)', 'dl(mm)']
        temp = '%-10s %-8s %-8s %-12s %-12s %-12s'
        cont.append(temp%tuple(key)); cont.append('-'*63)
        temp = '%-10s %-8s %-8s %-12.4f %-12.4f %-12.4f'
        for item in [i for i in info if i[0]=='distance']:
            cont.append(temp%tuple(item))

    if info and 'level' in [i[0] for i in info]:
        cont.append('\n水准观测及精度:')
        key = ['measure', 'from', 'to', 'mea(m)', 'value(m)', 'dl(mm)']
        temp = '%-10s %-8s %-8s %-12s %-12s %-12s'
        cont.append(temp%tuple(key)); cont.append('-'*63)
        temp = '%-10s %-8s %-8s %-12.4f %-12.4f %-12.4f'
        for item in [i for i in info if i[0]=='level']:
            cont.append(temp%tuple(item))
            
    if sta>0: cont.append('\n后验单位中误差:%.4f'%sta**0.5)
    return '\n'.join(cont)

# 绘制网点图
def plot(pts, info, sta):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
        
    dicpts = dict([(i[0], i) for i in pts])
    meas = []
    for mea in info:
        if mea[0]=='angle':
            _, p102, m, v, dl = mea
            p1, p0, p2 = p102.split('-')
            p01, p02 = sorted([p0,p1]), sorted([p0,p2])
            if not p01 in meas: meas.append(p01)
            if not p02 in meas: meas.append(p02)
        if mea[0]=='distance' or mea[0]=='direction':
            p12 = sorted(mea[1:3])
            if not p12 in meas: meas.append(p12)
    for p1, p2 in meas:
        p1, p2 = dicpts[p1], dicpts[p2]
        ax.plot([p1[2],p2[2]], [p1[1],p2[1]], 'gray')
        
    for name, x, y, _, is_true in [i[:5] for i in pts]:
        color = 'ro' if is_true else 'go'
        ax.plot([y], [x], color)
        ax.text(y, x, name, fontsize=12, color='black')
    
    ax.set_title('Points Plot')
    ax.set_xlabel('E')
    ax.set_ylabel('N')
    plt.show()
