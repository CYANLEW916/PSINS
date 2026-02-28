# 实验设置 (Experimental Setup)

本文档梳理仿真实验的完整设置流程，涵盖飞行程序选择、轨迹生成、传感器噪声建模与故障注入四个环节。所有参数均可在 `config.m` 中统一配置，仿真主流程由 `main_simulation.m` 驱动。

---

## 1. RNP AR 飞行程序的选择

### 1.1 程序背景

本仿真选用深圳宝安国际机场 **NLG RNP AR 进近程序** 作为飞行测试场景。RNP AR（Required Navigation Performance - Authorization Required）进近程序要求机载导航系统提供高精度的水平引导能力，是验证冗余惯导系统故障检测性能的典型工况。

该程序定义在 `main_simulation.m:177-186` 的 `generate_trajectory()` 函数中，同时在 `demos/test_SINS_trj_NLG` 中有独立的演示脚本。

### 1.2 航路点定义

进近程序包含 7 个航路点，从 IAF（起始进近定位点）到 RF_End（最终进近段终点）：

| 序号 | 航路点名称 | 角色 | 纬度 (°N) | 经度 (°E) | 高度 (ft MSL) | 说明 |
|------|-----------|------|-----------|-----------|--------------|------|
| 1 | NLG | IAF | 22.5317 | 113.5617 | 4900 | 起始进近定位点 |
| 2 | SZ928 | Fly-by | 22.4867 | 113.6517 | 2600 | 飞越航路点 |
| 3 | SZA34 | Fly-by | 22.4417 | 113.7400 | 2300 | 飞越航路点 |
| 4 | SZ924 | IF | 22.4683 | 113.7967 | 1300 | 中间进近定位点 |
| 5 | SZ923 | Fly-by | 22.4783 | 113.8183 | 1100 | 飞越航路点 |
| 6 | SZ921 | Fly-by | 22.4867 | 113.8350 | 1100 | 飞越航路点 |
| 7 | SZ920 | RF_End | 22.5333 | 113.8517 | 1100 | 最终进近段终点 |

> 源码位置：`main_simulation.m:179-185`，对应 `makeWP()` 函数调用。

### 1.3 航路点角色说明

- **IAF (Initial Approach Fix)**：进近程序的入口点，飞机从航路过渡到进近阶段。
- **IF (Intermediate Fix)**：中间进近定位点，衔接起始进近段与最终进近段。
- **Fly-by**：飞越式航路点，飞机在接近时提前转弯，不必精确飞越该点上空。
- **RF_End (Radius-to-Fix End)**：圆弧至定位点段的终点，标志最终进近段结束。

### 1.4 航段速度配置

每段航路的目标地速按进近阶段逐步递减（`main_simulation.m:172`）：

| 航段 | 起止航路点 | 目标速度 (kts) |
|------|-----------|---------------|
| Leg 1 | NLG → SZ928 | 210 |
| Leg 2 | SZ928 → SZA34 | 200 |
| Leg 3 | SZA34 → SZ924 | 190 |
| Leg 4 | SZ924 → SZ923 | 170 |
| Leg 5 | SZ923 → SZ921 | 150 |
| Leg 6 | SZ921 → SZ920 | 140 |

速度单位在代码中通过 `glv.kn`（1 knot ≈ 0.5144 m/s）转换为 m/s。

### 1.5 飞行动力学参数

| 参数 | 值 | 说明 |
|------|----|------|
| `accelTime` | 60 s | 从静止加速到首段目标速度的参考时间 |
| `decelTime` | 40 s | 减速至停止的参考时间 |
| `turnRate` | 3 °/s | 标准转弯速率 |
| `turnRollTime` | 6 s | 协调转弯的滚转过渡时间 |
| `pitchRate` | 0.5 °/s | 俯仰姿态变化速率 |

> 源码位置：`main_simulation.m:172-175`

---

## 2. 仿真轨迹的生成

### 2.1 生成流水线

轨迹生成遵循以下流水线，全部在 `main_simulation.m:169-275` 的 `generate_trajectory()` 函数中实现：

```
航路点定义 (W)
    │
    ▼
waypointLegs(W, legSpeeds)     ── 计算各航段的航向、俯仰、距离、时长
    │
    ▼
trjsegment()                   ── 逐段构建运动廓线 (加速/匀速/转弯/俯仰)
    │
    ▼
trjsimu(avp0, seg.wat, 0.01)   ── 以 0.01s 步长数值积分生成 IMU 增量数据
    │
    ▼
降采样至 cfg.ts = 0.02s         ── 角增量/速度增量累加，时间戳取末尾
    │
    ▼
trj.imu (Nx7)                  ── 最终仿真轨迹数据
```

### 2.2 航段几何计算 — `waypointLegs()`

`waypointLegs()` 函数（`main_simulation.m:296-318`）接收航路点序列和速度数组，计算相邻航路点间的飞行几何参数：

- **course**：航段航向角（弧度），由东向/北向分量的反正切计算：`atan2(ΔE, ΔN)`
- **pitch**：航段俯仰角（弧度），由高度差和水平距离计算：`atan2(Δalt, horizontal)`
- **horizontal**：水平距离 (m)，基于 PSINS 的 `pos2dxyz()` 函数进行大地坐标到 ENU 局部坐标的转换
- **vertical**：高度变化 (m)
- **time**：航段飞行时长 = `horizontal / speed`
- **speed**：航段目标速度 (m/s)

### 2.3 运动分段 — `trjsegment()`

`trjsegment()`（位于 `base/base1/trjsegment.m`）将飞行轨迹分解为基本运动单元，每个分段通过 WAT 矩阵（Waypoint And Trajectory）描述，每行包含 8 个参数：`[持续时间, 速度, ωx, ωy, ωz, ax, ay, az]`。

本仿真中使用的分段类型：

| 分段类型 | 用途 | 关键参数 |
|---------|------|---------|
| `init` | 初始化，设定初始速度 | 速度值 |
| `uniform` | 匀速直线飞行 | 持续时间 |
| `accelerate` | 恒加速度飞行 | 持续时间, 加速度 |
| `deaccelerate` | 恒减速度飞行 | 持续时间, 减速度 |
| `headup` | 抬头爬升 | 持续时间, 俯仰速率 |
| `headdown` | 低头下降 | 持续时间, 俯仰速率 |
| `coturnleft` | 协调左转弯（含滚转过渡） | 转弯时间, 转弯速率, 滚转时间 |
| `coturnright` | 协调右转弯（含滚转过渡） | 转弯时间, 转弯速率, 滚转时间 |

### 2.4 各航段的运动序列

对于每一段航路（第 2 段起），运动构建顺序为（`main_simulation.m:212-239`）：

1. **航向转弯**：若航向变化量 > 0.5°，插入协调转弯分段（`coturnleft` 或 `coturnright`）
2. **俯仰调整**：若俯仰角变化量 > 0.1°，插入爬升/下降分段（`headup` 或 `headdown`）
3. **速度调整**：若目标速度发生变化，插入加速/减速分段
4. **匀速飞行**：以目标速度匀速飞行至下一航路点

飞行结束后，恢复水平飞行并减速至零（`main_simulation.m:241-250`）。

### 2.5 轨迹仿真器 — `trjsimu()`

`trjsimu()`（位于 `base/base1/trjsimu.m`）根据 WAT 矩阵进行数值积分，逐步计算飞行器的姿态、速度和位置，并生成对应的 IMU 传感器增量数据。

**初始状态** `avp0` 的构造（`main_simulation.m:189`）：

```matlab
avp0 = [[0; 0; legs(1).course]; [0; 0; 0]; [W(1).lat; W(1).lon; W(1).alt]];
%        姿态 [俯仰=0, 滚转=0, 航向]    速度 [0,0,0]    位置 [纬度, 经度, 高度]
```

**核心积分过程**：每个时间步内，利用方向余弦矩阵 (DCM) 将轨迹帧运动转换为机体帧，生成陀螺角增量 `dθ` 和加速度计速度增量 `dv`。同时考虑地球自转 (`wnin`) 和重力加速度 (`gcc`) 的影响。

### 2.6 采样与降采样

- **生成采样率**：0.01s（100 Hz），由 `trjsimu(avp0, seg.wat, 0.01, 1)` 指定
- **目标采样率**：`cfg.ts = 0.02s`（50 Hz）
- **降采样方法**（`main_simulation.m:255-264`）：每 `ratio = ts/0.01 = 2` 个原始样本合并为一个：
  - 角增量和速度增量：**累加**（保持增量语义）
  - 时间戳：取最后一个样本的时间

### 2.7 输出数据格式

```
trj.imu  — N×7 矩阵: [dθ_x, dθ_y, dθ_z, dv_x, dv_y, dv_z, t]
trj.avp  — N×10 矩阵: [pitch, roll, yaw, vE, vN, vU, lat, lon, alt, t]
trj.avp0 — 初始状态向量
trj.ts   — 采样周期
```

### 2.8 轨迹缓存机制

为避免重复生成，轨迹数据采用文件缓存（`main_simulation.m:26-34`）：

```matlab
trjPath = fullfile(pwd, 'data', 'trj_NLG_approach.mat');
if exist(trjPath, 'file')
    trj = trjfile(trjPath);      % 直接加载已有轨迹
else
    trj = generate_trajectory(ts); % 首次运行时生成
    trjfile(trjPath, trj);        % 保存供后续使用
end
```

---

## 3. 传感器噪声的设置

### 3.1 冗余传感器架构

本仿真模拟 **2 套 INS + 1 套 ISIS** 共 9 个传感器通道的冗余配置：

| 子系统 | 通道索引 | 传感器编号 | 精度等级 |
|--------|---------|-----------|---------|
| INS1 | 1, 2, 3 | X, Y, Z 轴 | 高精度 |
| INS2 | 4, 5, 6 | X, Y, Z 轴 | 高精度 |
| ISIS | 7, 8, 9 | X, Y, Z 轴 | 低精度（备份） |

这种异构冗余配置是 SWGLT 算法的核心验证场景——ISIS 噪声比 INS 高一个数量级，传统 GLT 方法假设等精度传感器，无法有效处理这种差异。

### 3.2 噪声参数

传感器噪声标准差定义在 `config.m:17-23`，对应论文 Table 1：

**陀螺仪偏置不稳定性**：

| 传感器 | 原始值 | 单位转换 | σ (rad/s) |
|--------|--------|---------|-----------|
| INS1 / INS2 | 0.01 deg/h | × π/(180×3600) | 4.848×10⁻⁷ |
| ISIS | 0.1 deg/h | × π/(180×3600) | 4.848×10⁻⁶ |

**加速度计偏置不稳定性**：

| 传感器 | 原始值 | 单位转换 | σ (m/s²) |
|--------|--------|---------|----------|
| INS1 / INS2 | 1×10⁻⁴ g | × 9.8 | 9.8×10⁻⁴ |
| ISIS | 5×10⁻³ g | × 9.8 | 4.9×10⁻² |

> ISIS 的陀螺噪声是 INS 的 **10 倍**，加速度计噪声是 INS 的 **50 倍**。

### 3.3 噪声向量组装

噪声标准差组装为 1×9 向量，与 9 通道量测矩阵的列对齐（`config.m:26-32`）：

```matlab
cfg.sigma_gyro = [σ_ins, σ_ins, σ_ins,  σ_ins, σ_ins, σ_ins,  σ_isis, σ_isis, σ_isis]
%                 ├── INS1 (1-3) ──┤  ├── INS2 (4-6) ──┤  ├── ISIS (7-9) ──┤

cfg.sigma_acc  = [σ_ins, σ_ins, σ_ins,  σ_ins, σ_ins, σ_ins,  σ_isis, σ_isis, σ_isis]
```

### 3.4 配置矩阵 H

配置矩阵 `H`（`config.m:36`）描述 9 个传感器与 3 个物理量（X/Y/Z 轴）的量测关系：

```
        X  Y  Z
INS1-X [1  0  0]
INS1-Y [0  1  0]
INS1-Z [0  0  1]
INS2-X [1  0  0]     H = repmat(eye(3), 3, 1)   ——  9×3 矩阵
INS2-Y [0  1  0]
INS2-Z [0  0  1]
ISIS-X [1  0  0]
ISIS-Y [0  1  0]
ISIS-Z [0  0  1]
```

每个物理轴被 3 个传感器（INS1、INS2、ISIS）冗余测量，奇偶空间维度 = 9 - 3 = **6**。

### 3.5 量测方程

传感器量测模型为（`main_simulation.m:64-71`）：

```
Z = H × x_true + noise
```

其中：
- `x_true` 为 3 维真实物理量（角速率或比力），从轨迹数据提取
- `noise` 为各通道独立的零均值高斯白噪声：`noise_i ~ N(0, σ_i²)`

### 3.6 量测生成流程

完整的量测生成过程（`main_simulation.m:49-73`）：

```
1. 提取参考信号
   gyro_ref = imu_ref(:,1:3) / ts    % 角增量 → 角速率 (rad/s)
   acc_ref  = imu_ref(:,4:6) / ts    % 速度增量 → 比力 (m/s²)

2. 生成独立噪声（每个子系统各自独立）
   noise_ins1_gyro = σ_ins_gyro × randn(N, 3)
   noise_ins2_gyro = σ_ins_gyro × randn(N, 3)
   noise_isis_gyro = σ_isis_gyro × randn(N, 3)
   （加速度计同理）

3. 组装 9 通道量测矩阵
   Z_gyro = [gyro_ref + noise_ins1,  gyro_ref + noise_ins2,  gyro_ref + noise_isis]
   Z_acc  = [acc_ref + noise_ins1,   acc_ref + noise_ins2,   acc_ref + noise_isis]
   （均为 N×9 矩阵）
```

随机数种子固定为 `rng(42)`（`main_simulation.m:15`），确保结果可复现。

---

## 4. 故障注入的设定

### 4.1 故障条件总览

仿真设置 3 种故障条件，分别针对不同传感器、不同故障类型和不同量级，定义在 `config.m:41-59`：

| 条件 | 目标传感器 | 通道索引 | 量测类型 | 故障类型 | 时间区间 (s) | 故障幅值 |
|------|-----------|---------|---------|---------|-------------|---------|
| 1 | INS1 Y轴陀螺 | 2 | gyro | 硬故障（阶跃） | [120, 135]、[440, 455]、[595, 610] | 0.5、1.0、2.0 deg/h |
| 2 | INS2 X轴加速度计 | 4 | acc | 硬故障（阶跃） | [100, 115]、[420, 435]、[575, 590] | 0.005g、0.01g、0.02g |
| 3 | INS1 Z轴陀螺 | 3 | gyro | 软故障（斜坡） | [260, 310] | 0.02 deg/h/s |

### 4.2 传感器通道索引映射

`sensor_idx` 对应 9 通道量测矩阵 Z 的列索引：

| 索引 | 子系统 | 物理轴 | 故障条件 |
|------|--------|-------|---------|
| 1 | INS1 | X | — |
| 2 | INS1 | Y | **条件 1** |
| 3 | INS1 | Z | **条件 3** |
| 4 | INS2 | X | **条件 2** |
| 5 | INS2 | Y | — |
| 6 | INS2 | Z | — |
| 7 | ISIS | X | — |
| 8 | ISIS | Y | — |
| 9 | ISIS | Z | — |

### 4.3 硬故障（阶跃偏差）注入

硬故障在指定时间窗口内向传感器通道叠加**恒定偏差**（`inject_fault.m:20-29`）：

```
Z_faulty(t ∈ [t_start, t_end], idx) = Z(t, idx) + magnitude
```

每个故障条件可包含多个不连续的时间窗口，每个窗口有独立的幅值，模拟不同严重程度的传感器突变故障。

**条件 1 示例**（INS1 Y轴陀螺）：
- 第 1 段：t ∈ [120, 135]s，叠加 0.5 deg/h 偏差
- 第 2 段：t ∈ [440, 455]s，叠加 1.0 deg/h 偏差
- 第 3 段：t ∈ [595, 610]s，叠加 2.0 deg/h 偏差

**条件 2 示例**（INS2 X轴加速度计）：
- 第 1 段：t ∈ [100, 115]s，叠加 0.005g 偏差
- 第 2 段：t ∈ [420, 435]s，叠加 0.01g 偏差
- 第 3 段：t ∈ [575, 590]s，叠加 0.02g 偏差

### 4.4 软故障（斜坡偏差）注入

软故障在单一时间窗口内向传感器通道叠加**线性增长偏差**（`inject_fault.m:31-42`）：

```
Z_faulty(t ∈ [t_start, t_end], idx) = Z(t, idx) + rate × (t - t_ref)
```

**条件 3 示例**（INS1 Z轴陀螺）：
- 时间窗口：t ∈ [260, 310]s（持续 50 秒）
- 斜率：0.02 deg/h/s
- 参考时间：t_ref = 260s
- 故障增长：从 0 线性增长到 0.02 × 50 = 1.0 deg/h

软故障模拟传感器性能的渐进退化，对故障检测算法的灵敏度提出更高要求。

### 4.5 故障掩码 (fault_mask)

仿真主循环中根据故障条件生成布尔型故障掩码（`main_simulation.m:107-115`），用于后续性能评估：

- **硬故障**（条件 1、2）：合并所有时间窗口的掩码，`fault_mask = ∪(t ∈ [t_start_i, t_end_i])`
- **软故障**（条件 3）：单一窗口掩码，`fault_mask = (t ∈ [t_start, t_end])`

### 4.6 性能评估指标

每种故障条件下，三种算法（GLT、WGLT、SWGLT）的检测结果通过 `evaluate_performance.m` 计算以下指标：

| 指标 | 公式 | 说明 |
|------|------|------|
| FDR (故障检测率) | TP / n_fault | 故障期间正确检测的比例 |
| FAR (虚警率) | FP / n_normal | 正常期间错误报警的比例 |
| Accuracy (准确率) | (TP + TN) / N | 整体判别正确率 |
| Delay (检测延迟) | t_first_detect − t_fault_start | 仅软故障：从故障起始到首次检测的时间 |

其中：
- **TP**：故障期间检测函数超过阈值的样本数
- **FP**：正常期间检测函数超过阈值的样本数
- **TN**：正常期间检测函数未超过阈值的样本数
- 检测阈值：GLT/WGLT 使用 χ²(6, P_FA=0.01) ≈ 16.81；SWGLT 使用自适应阈值

---

## 附录：仿真执行流程总览

```
main_simulation.m 执行流程：

1. 初始化
   glvs;  rng(42);  cfg = config();

2. 生成/加载轨迹
   trj = generate_trajectory(ts)  或  trj = trjfile(trjPath)

3. 生成冗余传感器量测
   Z_gyro, Z_acc  ← 参考信号 + 独立高斯噪声 (N×9)

4. 遍历 3 种故障条件
   ├── 4a. 注入故障          → inject_fault(Z, t, fault_cfg, type)
   ├── 4b. 运行 GLT          → fdi_glt(Z_work, σ_uni, cfg)
   ├── 4c. 运行 WGLT         → fdi_wglt(Z_work, σ_vec, cfg)
   ├── 4d. 运行 SWGLT        → fdi_swglt(Z_work, σ_vec, cfg)
   ├── 4e. 评估性能          → evaluate_performance(FD, T_D, fault_mask, ...)
   └── 4f. 绘制结果          → plot_results(...)

5. 输出统计结果与可视化图表
```
