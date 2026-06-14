# Window 停止规则与步长诊断说明

## 为什么需要 window 停止规则

原始 ScaledGD 停止规则只看相邻两步 objective 的相对变化。诊断图显示，在若干实验中 objective 曲线会很早变平，但矩阵估计 `M_k` 仍在持续移动，`||M_k-M0||_F` 仍然下降。因此，只看单步 objective change 可能过早停止。

单步 Frobenius movement 也不稳健。比如默认 `p=10,q=10,m=200` 设置下，使用 `rel_frob < 1e-3` 的单步规则会在第 289 步停止，但此时 `||M-M0||_F = 1.819`，明显差于继续迭代后的结果。

## 两种 `scaledgd-stop-window` 解释

在非 `window` 规则下：

```text
scaledgd-stop-window = W
```

表示停止条件需要连续满足 `W` 次。例如 `stop-rule=frob, stop-window=10` 表示连续 10 次 Frobenius movement 都低于阈值才停。

在 `stop-rule=window` 下：

```text
scaledgd-stop-window = W
```

表示移动窗口长度为 `W`。程序会计算最近 `W` 步的平均 objective 下降率和平均相对 Frobenius movement。只有两个平均量都低于阈值时才停。

## 当前采用的 window 统计量

设当前迭代为第 `k` 步，窗口长度为 `W`：

```text
window_rel_objective =
  (obj_{k-W} - obj_k) / ((1 + abs(obj_{k-W})) * W)

window_rel_frob =
  mean(||M_j - M_{j-1}||_F / (1 + ||M_{j-1}||_F), j = k-W+1,...,k)
```

停止条件为：

```text
window_rel_objective < scaledgd_tol
window_rel_frob      < scaledgd_frob_tol
```

这个规则的优点是，它不依赖真实 `M0`，可以在实际 simulation 中使用。

## 关于窗口长度

增大 `scaledgd-stop-window` 本身不会改变迭代路径。它只会改变什么时候停止。

如果原规则停得太早，增大 window 可能改善最终估计，因为算法会多跑一段。但如果已经达到 `maxit`，继续增大 window 不会改善结果，应该提高 `maxit` 或适当调整阈值。

在默认 `m=200` 诊断中：

- `window=50/100, tol=1e-6, maxit=500` 没有提前停，因为 window objective decrease 仍明显高于阈值。
- `window=100, tol=3e-6, maxit=1000` 在第 996 步停止，Frobenius error 接近 2000 步参考结果。

## 关于步长

步长太小会显著拖慢收敛。在 `p=30,q=25,m=500,cov-scale=0.2` 实验中：

- `eta=0.1` 和 `eta=0.25` 在 5000 步内没有满足严格 window 停止。
- `eta=0.5` 在 3378 步收敛。
- `eta=1` 在 2060 步收敛，并得到本组实验中最小的 Frobenius error。
- `eta=1.5` 进一步减少迭代数到 1540，但最终误差略高于 `eta=1`。
- `eta=2` 和 `eta=4` 虽然更快触发停止，但 backtracking 明显增加，最终误差没有改善。

因此，在当前实验设置下，推荐从 `eta=1` 开始；若只关心更少迭代，可尝试 `eta=1.5`，但不建议盲目继续增大步长。

