import sys; sys.path.append('../')
import measurefit as mfit

print('\n水准网平差演示，武汉大学，误差理论于测量平差基础（第三版） 114页 例7-1', '\n'+'='*98)
points = [
    ('A', 0, 0, 237.483, True),
    ('B', 0, 0, 0, False),
    ('C', 0, 0, 0, False),
    ('D', 0, 0, 0, False)
]

# 高差，路线距离（负数，km为单位），或观测精度
measures = [
    ('lev', 'A', 'B', 5.835, -3.500),
    ('lev', 'B', 'C', 3.782, -2.700),
    ('lev', 'A', 'C', 9.640, -4.000),
    ('lev', 'D', 'C', 7.384, -3.000),
    ('lev', 'A', 'D', 2.270, -2.500)
]

pts, info, sta = mfit.solve(points, measures, db=1, dk=1, da=1)
print(mfit.report(pts, info, sta, fmt='txt'))


print('\n测向网平差演示，武汉大学，误差理论于测量平差基础（第三版） 120页 例7-3', '\n'+'='*98)
points = [
    ('A', 8986.68, 5705.03, 0, True),
    ('B', 13737.37, 10501.92, 0, True),
    ('C', 6642.27, 14711.75, 0, True),
    ('D', 10122.12, 10312.47, 0, False)
]

measures = [
    ('station', 'D', None),
    ('dir', 'D', 'C', 0, None),
    ('dir', 'D', 'A', 127.48412, None),
    ('dir', 'D', 'B', 234.39243, None),
    
    ('station', 'C', None),
    ('dir', 'C', 'A', 0, None),
    ('dir', 'C', 'D', 23.45162, None),
    
    ('station', 'A', None),
    ('dir', 'A', 'B', 0, None),
    ('dir', 'A', 'D', 30.52440, None),
    ('dir', 'A', 'C', 59.18490, None),
    
    ('station', 'B', None),
    ('dir', 'B', 'D', 0, None),
    ('dir', 'B', 'A', 42.16391, None)
]

pts, info, sta = mfit.solve(points, mfit.dms2rs(measures), db=1, dk=1, da=1)
print(mfit.report(pts, info, sta, fmt='txt'))


print('\n测角网平差演示，武汉大学，误差理论于测量平差基础（第三版） 123页 例7-4', '\n'+'='*98)
points = [
    ('A', 8986.68, 5705.03, 0, True),
    ('B', 13737.37, 10501.92, 0, True),
    ('C', 6642.27, 14711.75, 0, True),
    ('D', 10122.12, 10312.47, 0, False)
]

measures = [
    ('ang', 'A', 'D', 'B', 106.50422, None),
    ('ang', 'D', 'A', 'B', 30.52440, None),
    ('ang', 'D', 'B', 'A', 42.16391, None),
    ('ang', 'D', 'A', 'C', 28.26050, None),
    ('ang', 'A', 'D', 'C', 127.48412, None),
    ('ang', 'D', 'C', 'A', 23.45162, None)
]

pts, info, sta = mfit.solve(points, mfit.dms2rs(measures), db=1, dk=1, da=1)
print(mfit.report(pts, info, sta, fmt='txt'))


print('\n测角网平差演示，武汉大学，误差理论于测量平差基础（第三版） 140页 例7-9', '\n'+'='*98)
points = [
    ('A', 9684.28, 43836.82, 0, True),
    ('B', 10649.55, 31996.50, 0, True),
    ('C', 19063.66, 37818.86, 0, True),
    ('D', 17814.63, 49923.19, 0, True),
    ('P1', 13188.61, 37334.97, 0, False),
    ('P2', 15578.61, 44391.03, 0, False),
]

measures = [
    ('ang', 'A', 'P1', 'B', 126.14241, None),
    ('ang', 'P1', 'A', 'B', 23.39469, None),
    ('ang', 'P1', 'B', 'A', 30.05467, None),
    ('ang', 'A', 'P2', 'D', 117.22462, None),
    ('ang', 'P2', 'A', 'D', 31.26500, None),
    ('ang', 'P2', 'D', 'A', 31.10226, None),
    ('ang', 'P2', 'C', 'D', 22.02430, None),
    ('ang', 'C', 'P2', 'D', 130.03142, None),
    ('ang', 'P2', 'D', 'C', 27.53593, None),
    ('ang', 'P1', 'P2', 'A', 65.55008, None),
    ('ang', 'P1', 'A', 'P2', 67.02494, None),
    ('ang', 'P2', 'P1', 'A', 47.02114, None),
    ('ang', 'C', 'P2', 'P1', 46.38564, None),
    ('ang', 'C', 'P1', 'P2', 66.34547, None),
    ('ang', 'P1', 'C', 'P2', 66.46082, None),
    ('ang', 'P1', 'C', 'B', 29.58355, None),
    ('ang', 'C', 'P1', 'B', 120.08311, None),
    ('ang', 'C', 'B', 'P1', 29.52554, None)
]

pts, info, sta = mfit.solve(points, mfit.dms2rs(measures), da=1)
print(mfit.report(pts, info, sta, fmt='txt'))


print('\n测边网平差演示，武汉大学，误差理论于测量平差基础（第三版） 145页 例7-10', '\n'+'='*98)
points = [
    ('A', 53743.136, 61003.826, 0, True),
    ('B', 47943.002, 66225.854, 0, True),
    ('C', 40049.229, 53782.790, 0, True),
    ('D', 36924.728, 61027.086, 0, True),
    ('P1', 48580.270, 60500.505, 0, False),
    ('P2', 48681.390, 55018.279, 0, False),
    ('P3', 43767.223, 57968.593, 0, False),
    ('P4', 40843.219, 64867.875, 0, False)
]

measures = [
    ('dist', 'P1', 'B', 5760.706, None),
    ('dist', 'P1', 'A', 5187.342, None),
    ('dist', 'P2', 'A', 7838.880, None),
    ('dist', 'P2', 'P1', 5483.158, None),
    ('dist', 'P2', 'P3', 5731.788, None),
    ('dist', 'P2', 'C', 8720.162, None),
    ('dist', 'P3', 'C', 5598.570, None),
    ('dist', 'P3', 'D', 7494.881, None),
    ('dist', 'P3', 'P4', 7493.323, None),
    ('dist', 'P3', 'P1', 5438.382, None),
    ('dist', 'P4', 'D', 5487.073, None),
    ('dist', 'P4', 'P1', 8884.587, None),
    ('dist', 'P4', 'B', 7228.367, None)
]

pts, info, sta = mfit.solve(points, measures, maxiter=3, accu=0.1, db=3, dk=1, lw=1)
print(mfit.report(pts, info, sta, fmt='txt'))


print('\n带有距离，方向约束的测角网平差演示，武汉大学，误差理论于测量平差基础（第三版） 176页 例8-2', '\n'+'='*98)
points = [
    ('A', 2794005.704, 19433831.155, 0, True),
    ('B', 2802234.190, 19437826.220, 0, True),
    ('C', 2804773.909, 19432985.959, 0, False),
    ('D', 2805958.639, 19426570.796, 0, False),
    ('E', 2799571.971, 19430754.937, 0, False),
    ('F', 2798372.250, 19423925.543, 0, False),
    ('G', 2793886.720, 19428172.793, 0, False)
]

measures = [
    ('ang', 'B', 'A', 'E', 54.492957, None),
    ('ang', 'E', 'B', 'A', 43.281822, None),
    ('ang', 'B', 'E', 'A', 81.421162, None),
    ('ang', 'E', 'B', 'C', 48.190060, None),
    ('ang', 'B', 'C', 'E', 85.313721, None),
    ('ang', 'C', 'E', 'B', 46.092035, None),
    ('ang', 'D', 'C', 'E', 77.145904, None),
    ('ang', 'C', 'D', 'E', 46.182193, None),
    ('ang', 'C', 'E', 'D', 56.263883, None),
    ('ang', 'E', 'D', 'F', 52.271149, None),
    ('ang', 'D', 'F', 'E', 60.484761, None),
    ('ang', 'D', 'E', 'F', 66.435930, None),
    ('ang', 'E', 'F', 'G', 56.313621, None),
    ('ang', 'F', 'G', 'E', 67.514911, None),
    ('ang', 'G', 'E', 'F', 55.363570, None),
    ('ang', 'E', 'G', 'A', 64.220649, None),
    ('ang', 'G', 'A', 'E', 62.163629, None),
    ('ang', 'A', 'E', 'G', 53.211420, None),
    
    ('dir', 'B', 'E', 249.221017, 0),
    ('dist', 'D', 'C', 6523.643, 0)
]

pts, info, sta = mfit.solve(points, mfit.dms2rs(measures), maxiter=3, accu=0.1, db=3, dk=1, lw=1)
print(mfit.report(pts, info, sta, fmt='txt'))
mfit.plot(pts, info, sta)
