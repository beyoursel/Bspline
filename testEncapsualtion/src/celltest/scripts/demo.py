import numpy as np
from scipy.interpolate import BSpline, splprep, splev

# 假设你有一些控制点 P_ij 和其 x, y, z 坐标
control_points = np.array([
    [[0, 0, 0], [1, 0, 1], [2, 0, 2], [3, 0, 3]],
    [[0, 1, 1], [1, 1, 2], [2, 1, 3], [3, 1, 4]],
    [[0, 2, 2], [1, 2, 3], [2, 2, 4], [3, 2, 5]],
    [[0, 3, 3], [1, 3, 4], [2, 3, 5], [3, 3, 6]],
])

# 分别提取x, y, z坐标
x = control_points[:, :, 0]
y = control_points[:, :, 1]
z = control_points[:, :, 2]

# 使用 splprep 函数找到参数 u, v
tck_u, _ = splprep(x.T, s=0)
tck_v, _ = splprep(y.T, s=0)

# 指定需要插值的 (x, y)
x_new, y_new = 1.5, 1.5

# 找到对应的参数值 u, v
u_new = splev(x_new, tck_u)
v_new = splev(y_new, tck_v)

# 计算基函数值（这里简化，实际应根据 u_new 和 v_new 计算）
# 假设我们有基函数的计算函数
N_u = [BSpline.basis_element(tck_u[0][i:i+4], extrapolate=False)(u_new[0]) for i in range(len(tck_u[0])-4)]
M_v = [BSpline.basis_element(tck_v[0][i:i+4], extrapolate=False)(v_new[0]) for i in range(len(tck_v[0])-4)]

# 计算 z 值
z_new = sum(N_u[i] * sum(M_v[j] * z[i][j] for j in range(len(M_v))) for i in range(len(N_u)))

print(f'Interpolated z value at (x={x_new}, y={y_new}) is {z_new}')
