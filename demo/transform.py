import sys; sys.path.append('../')
import measurefit as mfit

meas = [
    ('station', 'TN1', 100),
    ('ts', 'TN1', 'TN2', 30.15319, 15.37591, 1000)
]

print('\n度分秒改弧度:')
meas = mfit.dms2rs(meas)
print(meas)

print('\n计算气象ppm:')
ppm = mfit.count_ppm(35, 0.5, 1000)
print(ppm)

print('\n加乘常数,ppm改正:')
meas = mfit.ts_adjusts(meas, b=1, k=1, ppm=ppm)
print(meas)

print('\n水平角，竖直角，斜距转水平角，高差，平距，并进行曲率改正，折光改正')
meas = mfit.ts2ahds(meas, db=1, dk=1, da=1, K=0.1)
print(meas)

print('\n平距投影到500m海拔')
meas = mfit.distprojs(meas, 500)
print(meas)

# 高效算法是基于axhd计算的，经典模型是基于dir, dist, lev计算的。
print('\n拆分成基础观测')
meas = mfit.splitmeas(meas)
print(meas)

print('\n一个完整的例子')

# 后方交汇
points = [
    ('TN2', 1000, 1000, 101, True),
    ('TN8', 1367.77286, 1117.43005, 102, True),
    ('TN9', 1294.76794, 1789.72043, 103, True),
    ('TN1', 893, 1603, 100, False)]

measures = [
    ('station', 'TN1', 100),
    ('ts', 'TN1', 'TN2', 0, 0.00163, 612.41862),
    ('ts', 'TN1', 'TN8', 0.598815, 0.00294, 678.73168),
    ('ts', 'TN1', 'TN9', 1.8317007, 0.00678, 442.39598)
]

print('\n预处理')
ppm = mfit.count_ppm(35, 0.5, 1000)
meas = mfit.ts_adjusts(measures, b=1, k=1, ppm=ppm)
meas = mfit.ts2ahds(meas, db=1, dk=1, da=1, K=0.1)
meas = mfit.distprojs(meas, 500)

print('\n高效算法到这里即可')
rst, sta = mfit.backward_x(points, meas)
print(mfit.report([rst], fmt='txt'))


print('\n经典算法要继续拆分成基础观测')
meas = mfit.splitmeas(meas)
pts, info, sta = mfit.solve(points, meas, db=1, dk=1, da=1)
print(mfit.report(pts, info, sta, fmt='txt'))



