from .base import *
from .util import *
from .io import *

def test():
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

    pts, info, sta = solve(points, dms2rs(measures), maxiter=3, accu=0.1, db=3, dk=1, lw=1)
    print(report(pts, info, sta))
    plot(pts, info, sta)