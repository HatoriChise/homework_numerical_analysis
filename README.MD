# 数值分析大作业 (Numerical Analysis Homework)

## 项目概述

本项目包含了数值分析课程的四次大作业，分别涉及不同数值方法的实现与应用。

## 项目结构

```
homework_numerical_analysis/
├── 2025数值计算大作业.pdf      # 原始作业要求文档
├── requirement.txt            # Python环境依赖
├── hw1/                     # 第一次作业：Romberg数值积分
│   ├── src/
│   └── picture/
├── hw2/                     # 第二次作业：三次样条插值 (C++)
│   ├── src/
│   ├── include/
│   ├── tests/
│   ├── CMakeLists.txt
│   └── picture/
├── hw3/                     # 第三次作业：Runge-Kutta方法求解微分方程
│   ├── src/
│   └── picture/
├── hw4/                     # 第四次作业：线性方程组求解器
│   ├── src/                 # 源文件目录
│   │   ├── solver/          # 求解器实现
│   │   └── *.cpp            # 主要源文件
│   ├── include/             # 头文件目录
│   │   └── solver/          # 求解器头文件
│   ├── build/               # 编译生成目录
│   ├── picture/             # 图片输出目录
│   └── CMakeLists.txt       # CMake构建配置
└── README.MD                # 本文件
```

## 各次作业说明

### HW1: Romberg数值积分算法

- **目的**: 实现Romberg积分算法，用于数值计算定积分
- **技术**: Python实现
- **主要内容**: 
  - 集成器类，实现Romberg积分方法
  - 收敛历史记录和可视化
  - 被积函数分析

### HW2: 三次样条插值

- **目的**: 实现三次样条插值算法
- **技术**: C++实现
- **主要内容**:
  - 自然边界条件的三次样条
  - 边界条件设置
  - 插值结果可视化

### HW3: Runge-Kutta方法求解常微分方程

- **目的**: 实现Runge-Kutta方法求解微分方程
- **技术**: Python实现
- **主要内容**:
  - 经典四阶Runge-Kutta算法实现
  - 与SciPy结果对比
  - 参数敏感性分析
  - 误差分析

### HW4: 线性方程组求解器

- **目的**: 实现线性方程组求解器，包含多种数值方法
- **技术**: C++实现，面向对象设计
- **主要内容**:
  - **自定义线性代数库**: 实现Matrix类用于矩阵操作
  - **LU分解求解器**: 直接法求解线性方程组Ax=b
  - **逐次超松弛(SOR)求解器**: 迭代法求解线性方程组
  - **抽象求解器基类**: 支持多种求解算法的统一接口
  - **系统方程类**: 管理系数矩阵、右端向量和求解过程
- **特点**:
  - 模块化设计，支持多种求解算法切换
  - 面向对象架构，代码结构清晰
  - 包含Python验证脚本用于结果对比
  - 支持任意大小的矩阵的LU分解和SOR迭代
  - 后续待补充基于稠密矩阵和稀疏矩阵的存储优化和求解器优化

## 环境配置

项目依赖Python和C++环境。

### Python环境配置

```bash
# 使用conda创建环境
conda create --name numerical_analysis --file requirement.txt
conda activate numerical_analysis
```

### C++环境配置

HW2和HW4使用C++，需要C++编译器和CMake。

```bash
# 在hw2目录下构建
cd hw2
mkdir build
cd build
cmake ..
make
```

## 运行各次作业

### HW1
```bash
cd hw1/src
python main.py
```

### HW2
```bash
cd hw2/build
make
./main
```

### HW3
```bash
cd hw3/src
python main.py
```

### HW4
```bash
cd hw4
mkdir build
cd build
cmake ..
make
./system_equation
```

## 依赖库

Python部分依赖以下库:
- numpy
- matplotlib
- scipy

C++部分使用标准库实现。

构建工具:
- CMake (用于HW2和HW4)
- C++17编译器

## 项目特色

1. **多语言实现**: 涵盖Python和C++两种编程语言
2. **理论与实践结合**: 算法实现与数值验证相结合
3. **可视化分析**: 提供收敛性、误差分析等可视化结果
4. **模块化设计**: 代码结构清晰，便于理解与扩展

## 作者

项目完成者: 谢振辉

## 日期

2025年11月11日数值分析课程

## 备注

HW4的Python验证脚本（true_result.py）位于`hw4/src/true_result.py`，使用NumPy和SciPy验证C++实现的正确性。项目完整实现了数值分析课程的全部四个作业，涵盖了数值积分、插值、微分方程求解和线性方程组求解等核心数值方法。