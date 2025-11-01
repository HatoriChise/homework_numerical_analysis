import mpmath as mp

mp.mp.dps = 25  # 设置十进制精度（25位足够）


def f(x):
    t = mp.pi * x / 50
    return (
        (4 * mp.pi**2 / 5)
        * (1 - mp.cos(t))
        * mp.fabs(mp.sin(t))
        * mp.sqrt(25 + 256 * mp.cos(t) ** 2)
    )


# 在 x=50 处分段（因为 |sin| 在此处不可导）
res = mp.quad(f, [0, 50, 100])
print(res)
