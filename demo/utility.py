import sys; sys.path.append('../')
import measurefit as mfit


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

print('\n后方交汇(经典算法)', '\n'+'='*98)
measx = mfit.ts2ahds(measures, 1, 1, 1, K=1)
measx = mfit.splitmeas(measx)
pts, info, sta = mfit.solve(points, measx, db=1, dk=1, da=1)
print(mfit.report(pts, info, sta, fmt='txt'))


print('\n后方交汇(高效算法)')
measx = mfit.ts2ahds(measures, 1, 1, 1, K=1) # 全站仪数据转水平角，高差，平距
rst, sta = mfit.backward_x(points, measx)
print(mfit.report([rst], fmt='txt'))




# 极坐标
points = [
    ('P1', 0, 0, 0, True),
    ('P2', 1000, 0, 0, True),
    ('Px', 1000, 577, 100, False)
]

# A站观测
meas_a = [
    ('station', 'P1', None),
    ('ts', 'P1', 'P2', 0, 0, 1000),
    ('ts', 'P1', 'Px', 0.52359, 0.0864, 1158.84),
    # 这里可以有更多未知点观测记录
]

# B站观测
meas_b = [
    ('station', 'P2', None),
    ('ts', 'P2', 'P1', 3.141592, 0, 1000),
    ('ts', 'P2', 'Px', 1.57079, 0.17151, 585.94),
    # 这里可以有更多未知点观测记录
]

print('\n前方交汇(经典算法)', '\n'+'='*98)
measx = mfit.ts2ahds(meas_a+meas_b, 1, 1, 1, K=1)
measx = mfit.splitmeas(measx)
pts, info, sta = mfit.solve(points, measx, db=1, dk=1, da=1)
print(mfit.report(pts, info, sta, fmt='txt'))


# 改正，转高差，平距
meas = mfit.ts_adjusts(meas_a, b=0, k=1, ppm=0)
meas = mfit.ts2ahds(meas, 1, 1, 1, K=1)
xyh_a = mfit.polar(points, meas)
print('\nA站极坐标(高效算法)')
print(mfit.report(xyh_a, fmt='txt'))


# 改正，转高差，平距
meas = mfit.ts_adjusts(meas_b, b=0, k=1, ppm=0)
meas = mfit.ts2ahds(meas, 1, 1, 1, K=1)
xyh_b = mfit.polar(points, meas)
print('\nB站极坐标(高效算法)')
print(mfit.report(xyh_b, fmt='txt'))

xyh_ab = mfit.merge_xyhs([xyh_a, xyh_b])
print('\nAB站极坐标成果融合（高效算法）')
print(mfit.report(xyh_ab, fmt='txt'))

print('\n')
meas = mfit.ts_adjusts(meas_a+meas_b, b=0, k=1, ppm=0) # 加乘常数改正，气象改正
meas = mfit.ts2ahds(meas, 1, 1, 1, K=1) # 全站仪数据转水平角，高差，平距

xyh_abx = mfit.forward_x(points, meas)
print('\nAB站前方交汇（高效算法，利用AB极坐标的结果进行融合）')
print(mfit.report(xyh_abx, fmt='txt'))

