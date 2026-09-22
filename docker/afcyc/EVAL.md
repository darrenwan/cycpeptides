# AfCyc Dockerfile 评估：科学 × 简练

**对象**：`docker/afcyc/Dockerfile`（镜像 `pdl1-afcyc`）  
**日期**：2026-09-22（修订后）  
**对照**：战役规格 S2c；环化在 `scripts/pdl1_afcyc_predict.py`。

## 总评（修订后）

| 维度 | 评分 | 一句话 |
|------|------|--------|
| **科学性** | **A−** | 单来源 ColabDesign pin + jax/haiku/numpy 钉死；权重外挂 |
| **简练性** | **B+** | 单次 JAX 安装；无 git fallback；构建入口已补齐 |

## 已落地的修改（相对初版 EVAL）

| 建议 | 状态 |
|------|------|
| P0 钉死 ColabDesign / 去掉双路径 | ✅ 仅 `COPY third_party/ColabDesign`；`fetch_colabdesign.sh` @ **v1.1.3**（`v1.1.1` 不存在） |
| P0 `build_academic_images.sh` | ✅ |
| P1 numpy / biopython 钉版本 | ✅ `1.26.4` / `1.84` |
| P1 SOURCES 记录 pin | ✅ `models/SOURCES.md` |
| P2 单次 JAX；缩短注释；去无用 WORKDIR | ✅；ColabDesign 用 `--no-deps` 避免先装 jax≥0.5 |
| 额外 | ✅ `.dockerignore`；`UV_HTTP_TIMEOUT=7200`；bootstrap 调用 fetch |

## 构建

```bash
bash scripts/build_academic_images.sh afcyc
# 等价：fetch_colabdesign → docker build -f docker/afcyc/Dockerfile -t pdl1-afcyc:latest .
```

## 仍注意

- 重建镜像后应用 `AFCYC_LIMIT=2 bash scripts/pdl1_smoke_afcyc.sh` 回归。  
- 主机 NVIDIA 驱动需支持 CUDA 12.1 runtime + jax 0.4.30。  
- AfCycDesign ≠ 官方 notebook 逐行拷贝；cyclic offset 在战役脚本中实现（`offset_type=2`）。
