## 数据格式介绍
measurefit 是完全面向过程的，没有任何自定义数据结构，类，全部采用python的基础容器实现。非常简单，高效。可以和json，yaml进行导入导出。

### 点数据
```python
points = [
    ('A', 8986.68, 5705.03, 0, True),
    ('B', 13737.37, 10501.92, 0, True),
    ('C', 6642.27, 14711.75, 0, True),
    ('D', 10122.12, 10312.47, 0, False)
]
```
点数据用一个python的tuple或list，分别是 Name, X, Y, H, Type。Type表示是否是已知点。


### 观测数据
```python
measures = [
    ('ang', 'B', 'A', 'E', 1.58635, None),  
    ('station', 'B', None),  
    ('dir', 'B', 'E', 0.221017, None),
    ('dist', 'D', 'C', 6523.643, None),
    ('lev', 'A', 'B', 12.75, None)
]
```

观测数据分多种，同样用一个tuple或list表示。
* **角度**：**[ang, p1, p2, p3, value, sta]**
* **测站**：**[station, name, h0]** 测站无需观测值，但携带了一个自由定向角，所以在方向观测之前，必须插入一个station记录。h0 是测站高程（概略即可），主要是用于做投影高程改正，如无需改正可给None。
* **方向**：**[dir, p1, p2, value, sta]**
* **平距**：**[dist, p1, p2, value, sta]**
* **水准**：**[lev, p1, p2, value, sta]**

**备注**：
* **单位**：上述距离观测值 m 为单位，角度观测弧度为单位。测量上惯用小数点表示度分秒，我们提供了dms2rs对观测数据进行转换。

* **精度**：sta表示观测精度，中误差。距离观测 mm 为单位，角度观测 ″ 为单位，如果指定则用给定精度进行平差，填None，求解器会利用原始测角测距精度进行推导。水准观测有时也会利用距离定权，用距离定权的时候，sta请输入距离的负数，m为单位，系统会用1km作为单位权。

* 观测精度给0表示约束条件，求解器将使用拉格朗日乘数进行约束平差求解。
* 
### 扩展观测数据
```python
measures = [
    ('station', 'TN1', None),
    ('ts', 'TN1', 'TN2', 0, 0.00163, 612.41862),
    ('axhd', 'TN1', 'TN2', 0.79, 0.16, 788,90, None, None, None)
]
```
扩展观测数据不是独立直接观测量，不能直接纳入求解器。是为了简化输入，需要使用转换函数进行转换。

* **全站仪观测**：**[ts, p1, p2, ax, ay, l]** 
* **方向高差平距观测**：**[axhd, ax, h, d dax, dh, d]** 

**转换方法**：  
* **ts2ahds** 函数将 **ts** 观测转换为 **axhd** 观测，这个过程中观测精度将自动计算并传播给 **axhd** 观测。
* **splitmeas** 将 **axhd** 观测分解为 **dir**, **lev**, **dist** 观测，进而求解。

### 解算点位成果
```python
rst = [('Px', 999.8374, 577.3316, 100.0019, False, 3.06, -0.65, -0.65, 1.45, 6.29)]
```
**成果点位数据**：**[Name, X, Y, H, Type, Dxx, Dxy, Dyx, Dyy, Dhh]** 后面分别是平面协方差以及高差方差。携带协方差信息，可以方便的计算点位精度，长短轴，以及多站同点成果融合。


## 通用求解函数
通用模型求解函数，传入测点数据，观测数据，即可进行平差。这部分的函数，在 [经典案例](classical.md) 中有详细的用法演示。

**`solve(pts, meas, maxiter=3, accu=0.1, db=1, dk=1, da=1, iter=False, lw=2)`**


#### 参数说明

- **pts**: 测点数据 参考 [数据格式](ioformat.md)
- **meas**: 观测数据，已知条件约束的精度给0，参考 [数据格式](doc/ioformat.md)
- **maxiter**: 迭代次数（类型: 整数，默认值: 3），规定次数内不收敛则抛出异常
- **accu**: 迭代精度控制（类型: 浮点数，默认值: 0.1，单位: mm）
- **db**: 测距固定误差（类型: 浮点数，默认值: 1，单位: mm）
- **dk**: 测距比例误差（类型: 浮点数，默认值: 1，单位: ppm）
- **da**: 测角误差（类型: 浮点数，默认值: 1，单位: ″）
- **iter**: 是否使用序贯方式求解（类型: bool，默认值: False）
- **lw**: 距离定权方式（类型: 整数, 1. `(db + dk * l) ^ 2`, 2. `db ^ 2 + (dk * l) ^ 2`， 默认值:2）

#### 返回值
- **pts**: 含有点位信息及点位精度信息的 list， 参考 [数据格式](doc/ioformat.md)
- **info**: 含有各个原始观测，观测平差值，观测精度的 list
- **sta**: 类型：float，单位权方差，开方后得单位权中误差

#### 异常

- **Exception**: 如果单步改正量大于 `accu`，且达到 `maxiter` 次数，则抛出异常。

## 高效求解函数
这部分的函数，在 [高效求解](utility.md) 中有详细的用法演示，此处只做简要说明。

**`backward_x(pts, meas, db=1, dk=1e-3, da=1, accu=1e-1, maxiter=3)`** 后方交汇，meas必须是从某个未知点看向若干已知点的观测记录。参数参考 solve

**`polar(pts, meas)`** 极坐标求解，观测必须是从已知点看向若干未知点，第一个是定向点

**`merge_xyh(xyhs)`** 融合同一个点的多个极坐标解算结果

**`merge_xyhs(xyhss)`** 融合多组极坐标解算结果，每一组可以包含若干点，不要求每组测点相同

**`forward_x(pts, meas)`** 前方交汇，meas必须是符合 polar 要求的观测，或几个符合polar要求观测的拼接

## 预处理函数
这部分的函数，在 [数据预处理](transform.md) 中有详细的用法演示，此处只做简要说明。

**`dms2r(v)`** 将小数形式的度分秒，转换成弧度值

**`dms2rs(meas)`** 传入观测列表，对里面的ts， angle， direction观测进行批量dms2r转换。

**`count_ppm(t, h, p)`** 根据温度，湿度，气压，计算距离改正比，单位ppm

**`ts_adjust(mea, b=0, k=1, ppm=0)`** 将一个 ts 观测进行加乘常数改正，以及气象ppm改正

**`ts_adjusts(meas, b=0, k=1, ppm=0)`** 将观测列表进行批量加乘常数，以及气象ppm改正

**`ts2ahd(mea, db, dk, da, K=0)`** 将一个 ts (ax, ay, l) 观测，进行球气差改正，转axhd (ax, h, d), 并进行精度传导

**`ts2ahds(meas, db, dk, da, K=0)`** 将观测列表中的 ts (ax, ay, l) 进行批量球气差改正，转axhd (ax, h, d), 并进行精度传导

**`def distproj(mea, h, h0)`** 将平距投影到统一水准面上，h 是测站测点平距高程，h0是目标投影高程

**`def distprojs(meas, h0)`** 将观测列表中的axdh （ax, h, d）的平距批量投影到 h0 高度。当前高度会利用meas中的station高程自动计算。

**`splitmeas`** 将axhd观测，分解成 dir, level, dist 观测，并分配精度。

## 输出函数
**`mapping(pts, info=None, sta=0.0)`** 将 solve 成果添加key，转成json

**`report(pts, info=None, sta=0.0)`** 将 solve 成果编制成可阅读的文本

**`plot(pts, info, sta)`** 将 solve 成果绘制出网点图