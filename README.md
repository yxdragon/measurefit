# MeasureFit
MeasureFit 是一个大地测量数据处理库，包括数据预处理，数据平差等功能。

### 特色
* 模型有很强的鲁棒性，可以同时输入多种观测，或约束条件。
* Python 语言实现，具有强大而灵活的io能力，成果渲染能力。
* 提供了一种多站极坐标融合成前方交会的高效，便捷算法。
* 没有使用任何数据结构，面向对象概念，全部基于python基础数据类型实现。

（暂无概略坐标计算功能，因为光电测量已经普及，多数情况下是有概略坐标的）

## 功能模块

### 平差模块
测量数据平差是大地测量中的核心任务之一，MeasureFit 提供了一个通用，鲁棒的平差模型：

* **支持多种观测**：水准观测，夹角观测，方向观测，距离观测，并支持多种观测的任意组合。
* **精度分析**：能够输出平差坐标，及精度分析，给出点位精度，以及各个观测量的平差精度。
* **成果输出**：能够输出json成果，txt报表，绘制图形。
* **支持条件约束**：支持含有已知高程，夹角，方向，距离的条件平差。
* **支持序贯求解**：支持观测量逐一纳入，采用卡尔曼滤波方式求解。

### 数据预处理
数据预处理是确保测量数据质量的关键步骤，MeasureFit 提供了以下功能来进行数据预处理：

* **加乘常数，气象改正**：可以对观测数据进行加乘常数改正，气象ppm改正。
* **球气差改正**：斜距转换高差平距的时候，可以进行球气差改正。
* **平距投影高程改正**：可以将平距投影到给定水准面上，做投影改正。

## 文档及示例

### [API文档](doc/api.md)
- **数据结构说明**：详细描述了数据结构的设计和使用方法。
- **函数列表**：列出了多数可用的函数及其功能说明。
- **功能说明**：对每个函数的功能进行了详细的解释和示例。

### [经典算例](doc/classical.md)
- **侧边**：侧边测量案例及其平差处理。
- **测角**：测角测量案例及其平差处理。
- **侧向**：侧向测量案例及其平差处理。
- **约束平差**：约束条件下的平差案例。

### [特色算例](doc/utility.md)
- **后方交汇**：后方交汇测量的高效实现。
- **前方交汇**：前方交汇测量的高效实现。
- **极坐标**：极坐标测量的高效实现。
- **极坐标融合**：极坐标融合测量的高效实现。

### [数据预处理](doc/transform.md)
- **加乘常数改正**：加乘常数改正的原理及实现。
- **气象ppm改正**：气象ppm改正的原理及实现。
- **球气差改正**：球气差改正的原理及实现。
- **高程面投影**：高程面投影的原理及实现。

## 算例演示
### 带有距离，方向约束的测角网平差演示，武汉大学，误差理论于测量平差基础（第三版） 176页 例8-2
```python
import measurefit as mfit

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
    # sta 给0代表约束，求解器会使用拉格朗日乘数进行约束求解
    ('dir', 'B', 'E', 249.221017, 0),
    ('dist', 'D', 'C', 6523.643, 0)
]

pts, info, sta = mfit.solve(points, mfit.dms2rs(measures), maxiter=3, accu=0.1, db=3, dk=1, lw=1)
print(mfit.report(pts, info, sta, fmt='txt'))
```

![netplot](https://github.com/user-attachments/assets/22ed6d03-71c1-4a4d-9fe8-c002efb7ec0c)

#### 点位信息及精度

| name | x            | y            | h     | type  | dxx       | dxy       | dyy       | dhh   |
|------|--------------|--------------|-------|-------|-----------|-----------|-----------|-------|
| A    | 2794005.7040 | 19433831.1550| 0.0000| True  | nan       | nan       | nan       | nan   |
| B    | 2802234.1900 | 19437826.2200| 0.0000| True  | nan       | nan       | nan       | nan   |
| C    | 2804773.9090 | 19432985.9625| 0.0000| False | 228.7628  | -35.5405  | 803.2532  | nan   |
| D    | 2805958.6406 | 19426570.7985| 0.0000| False | 865.9464  | 614.1075  | 1021.4715 | 0.0000|
| E    | 2799571.9725 | 19430754.9399| 0.0000| False | 42.2697   | 112.2752  | 298.2211  | 0.0000|
| F    | 2798372.2518 | 19423925.5449| 0.0000| False | 1084.6409 | 126.6167  | 995.1667  | 0.0000|
| G    | 2793886.7210 | 19428172.7942| 0.0000| False | 553.7872  | 47.2277   | 504.4942  | 0.0000|

#### 角度观测及精度

| measure | from_to | mea(°) | value(°) | dl(″) |
|---------|---------|--------|----------|-------|
| angle   | A-E-G   | 0.9312 | 53.3541  | 0.5447|
| angle   | B-A-E   | 0.9569 | 54.8247  | 0.3507|
| angle   | B-C-E   | 1.4927 | 85.5273  | 0.5777|
| angle   | B-E-A   | 1.4260 | 81.7031  | 0.3507|
| angle   | C-D-E   | 0.8082 | 46.3064  | 0.3116|
| angle   | C-E-B   | 0.8056 | 46.1558  | 0.5150|
| angle   | C-E-D   | 0.9851 | 56.4440  | 0.2767|
| angle   | D-C-E   | 1.3483 | 77.2496  | 0.5787|
| angle   | D-E-F   | 1.1647 | 66.7332  | 0.5503|
| angle   | D-F-E   | 1.0614 | 60.8136  | 0.4285|
| angle   | E-B-A   | 0.7587 | 43.4722  | 0.0000|
| angle   | E-B-C   | 0.8433 | 48.3168  | 0.4607|
| angle   | E-D-F   | 0.9155 | 52.4531  | 0.4459|
| angle   | E-F-G   | 0.9866 | 56.5265  | 0.5113|
| angle   | E-G-A   | 1.1234 | 64.3686  | 0.5473|
| angle   | F-G-E   | 1.1844 | 67.8638  | 0.5044|
| angle   | G-A-E   | 1.0869 | 62.2773  | 0.4951|
| angle   | G-E-F   | 0.9706 | 55.6098  | 0.5506|

#### 方向观测及精度

| measure  | from | to   | mea(°) | value(°) | dl(″) |
|----------|------|------|--------|----------|-------|
| direction| B    | E    | 249.3695| 249.3695 | 0.0000|

#### 边长观测及精度

| measure | from | to   | mea(m)   | value(m)  | dl(mm) |
|---------|------|------|----------|-----------|--------|
| distance| D    | C    | 6523.6430| 6523.6430 | 0.0000|

#### 后验单位中误差

后验单位中误差: 1.2033

