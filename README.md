# Support
如果觉得这个项目有用或对你有启发的话请给一颗⭐️，不胜感激。

If you find this project useful, please leave a ⭐️ for me. Thank you very much.

此为2026微生物应用创新设计大赛参赛作品湿实验验证部分
# 酵母木犀草素合成的 FBA / ecFBA 分析

这是我在项目中整理的一套代谢模型分析代码，主要想用 FBA 看木犀草素合成通路接入酵母后，会不会影响细胞生长，以及不同基因改造方案之间有什么区别。

项目以 yeastGEM 为基础，除了普通 FBA，我还尝试加入了一个比较简化的酶约束层，用来估计异源蛋白表达和 P450 酶可能带来的负担。这个模型主要用于前期比较和提出假设，不能直接代替真实实验结果。

## 1. 我做了哪些分析

- 比较野生型与加入木犀草素通路后的生长情况；
- 改变木犀草素最低输出通量，观察生长负担；
- 改变异源蛋白池和 P450 蛋白池，进行轻量 ecFBA 扫描；
- 按模块比较通路中不同部分的负担；
- 比较 `ARO7`、`PHA2` 敲除以及 `PheA^fbr` 引入后的结果；
- 画出通量需求与生长、不同模块负担的示意图。

## 2. 项目结构

```text
FBA/
├─ README.md                 # 项目说明和运行方法
├─ requirements.txt          # Python 依赖
├─ data/                     # 酶参数、模块设计和实验数据模板
├─ models/
│  └─ yeastGEM.xml           # 使用的酵母代谢模型
├─ src/                      # 建模、扫描和绘图代码
├─ results/                  # 运行后自动生成的 CSV 和图片
└─ docs/
   ├─ reports/               # 分析报告和阶段性总结
   └─ presentation/          # 汇报用 PPT
```

其中，`data/`、`models/` 和 `src/` 是复现分析所需的主要内容。`results/` 默认保持为空，运行脚本后会自动生成结果，所以不需要手动在里面新建文件。

## 3. 环境配置

建议使用 Python 3.10 及以上版本。在项目根目录打开 PowerShell，先用 `python --version` 检查版本；如果显示的是 Python 2.7，就不要继续使用这个命令，需要改用电脑中 Python 3 对应的启动命令或 `python.exe` 路径。

确认 Python 3 可以使用后，依次运行：

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
```

如果电脑上没有 Python 3.11，可以把第一行改成自己已经安装的 Python 3 版本，例如 `py -3.12 -m venv .venv`。如果 `py` 找不到已经安装的版本，也可以直接用 Python 3 的完整路径运行 `-m venv .venv`。

## 4. 运行方法

第一次运行时，可以先做一个普通 FBA，确认模型和环境能够正常读取：

```powershell
python src/run_fba.py
```

完整分析建议按下面的顺序运行：

```powershell
python src/fba_luteolin_burden.py
python src/run_fba.py
python src/run_ecfba_scan.py
python src/run_module_scan.py
python src/run_aro7_pha2_comparison.py
python src/run_phea_dual_function_analysis.py
python src/analyze_results.py
```

各脚本的作用如下：

| 脚本 | 作用 |
|---|---|
| `fba_luteolin_burden.py` | 完整通路的基础负担和产能场景 |
| `run_fba.py` | 普通 FBA 基线比较 |
| `run_ecfba_scan.py` | 扫描蛋白池、P450 池和木犀草素需求 |
| `run_module_scan.py` | 按通路模块进行比较 |
| `run_aro7_pha2_comparison.py` | 比较 `ARO7`、`PHA2` 等改造方案 |
| `run_phea_dual_function_analysis.py` | 分析双功能 `PheA^fbr` 模型 |
| `analyze_results.py` | 汇总部分结果并绘图 |

## 5. 输入和输出

主要输入文件：

- `models/yeastGEM.xml`：基础代谢网络；
- `data/enzymes.csv`：异源酶的分子量、`kcat` 和相对表达量等参数；
- `data/module_design.csv`：模块扫描时使用的反应组合；
- `data/experimental_results.csv`：实验数据模板。

运行后，结果会写入 `results/`，常见文件包括：

- `fba_luteolin_burden_results.csv`
- `fba_baseline.csv`
- `ecfba_scan.csv`
- `module_scan.csv`
- `aro7_pha2_comparison.csv`
- `phea_dual_function_analysis.csv`
- `figures/growth_vs_luteolin_demand.png`
- `figures/module_burden_comparison.png`

如果在 `data/experimental_results.csv` 中填入真实数据，再运行 `analyze_results.py`，还会生成 `model_experiment_gap.csv`，用于简单比较模型预测和实验结果。

## 6. 参数说明和默认值来源

为了避免把“为了扫描方便而设置的数值”误写成“实验测得的参数”，这里把参数来源分成四类：

- **模型内置值**：直接来自 `models/yeastGEM.xml`；
- **项目设计值**：由本项目采用的通路和基因改造方案决定；
- **扫描假设值**：为了比较趋势而人为设置，不是实测值；
- **文献依据值**：文献可以支持建模逻辑，但不一定能直接作为本模型的定量参数。

### 6.1 基础 FBA 和培养条件

| 参数 | 默认值 | 用来做什么 | 默认值来源 |
|---|---|---|---|
| `--model` | `models/yeastGEM.xml` | 指定基础酵母代谢模型 | 模型文件自身标记为 yeastGEM v9.0.2，对应 yeast-GEM 的 v9.0.2 发布版本 |
| `--growth-reaction` | 自动识别，当前为 `r_2111` | 指定最大化的生长反应 | `r_2111` 是本地 SBML 模型中已经设置好的目标反应，不是本项目新定义的参数 |
| `--glucose-exchange` | `r_1714` | 指定葡萄糖交换反应 | 来自 yeastGEM 的反应编号，模型中名称为 `D-glucose exchange` |
| `--oxygen-exchange` | `r_1992` | 指定氧气交换反应 | 来自 yeastGEM 的反应编号，模型中名称为 `oxygen exchange` |
| `--glucose-uptake` | `10.0 mmol/gDW/h` | 限制最大葡萄糖摄取速率 | 项目中的标准化培养条件假设，用于让不同场景在相同碳源供给下比较；不是本项目菌株的实测摄糖速率 |
| `--oxygen-uptake` | `1000.0 mmol/gDW/h` | 设置氧气摄取上限 | 数值足够大，相当于把氧气设为“不限量”；它表示有氧条件下不让氧气成为主要限制，不代表真实耗氧速率 |
| `--fns-type` | `p450` | 选择 FNS 反应类型；`p450` 消耗 NADPH 和 O₂，`2og` 使用 2-氧戊二酸 | 当前项目通路采用 FNSII/P450 路线，因此默认选 `p450`；`2og` 只作为替代路线对照 |
| `--include-side-reaction` | 关闭 | 是否加入 F3'H 将柚皮素转成圣草酚的旁路 | 默认只研究木犀草素主路线，所以关闭；开启后用于观察旁路竞争 |
| `--knockout` | `none` | 选择不敲除、敲除 `PHA2`、敲除 `ARO7` 或同时敲除 | 默认 `none` 用作未敲除基线；其他选项对应项目中要比较的基因改造方案 |
| `--luteolin-demands` | `fba_luteolin_burden.py`：`0,0.001,0.01,0.05`；`run_fba.py` 和 ecFBA 扫描另加入 `0.005` | 给木犀草素交换反应设置最低输出通量，比较产物要求升高后生长如何变化 | 扫描假设值，按“无产物要求—低—中—较高要求”设置，不是实验产率或滴度 |
| `--capacity-growth-floors` | `0.1,0.3` | 要求细胞至少保留 10% 或 30% 的 WT 生长，再求理论最大木犀草素通量 | 用于寻找模型能力边界的扫描假设，不是生物学阈值 |
| `--expression-atp-cost` | `0.0` | 给异源蛋白表达增加固定 ATP 水解负担 | 因为目前没有实测表达能耗，基础脚本默认不额外扣除 ATP，避免把猜测当成数据 |
| `--apply-capacity-modifiers` | 基础脚本默认关闭；模块表中的 M1–M4 和 M6 会启用 | 对部分中心代谢反应上限做倍数调整 | 这些倍数来自早期项目方案，只适合情景比较；yeastGEM 中很多原始上限是 `1000`，因此不能把缩放后的值解释成真实酶容量 |
| `--output` | 对应的 `results/*.csv` | 指定结果保存位置 | 项目文件组织约定，只改变输出文件，不改变模型结果 |

这里的通量单位通常写作 `mmol/gDW/h`，表示每克细胞干重、每小时的物质的量通量。需要特别注意：`luteolin demand = 0.001` 是最低通量约束，不是 `mg/L` 滴度，也不是发酵终点浓度。

### 6.2 轻量 ecFBA 参数

本项目没有直接建立完整 GECKO 模型，而是参考酶约束模型的基本思想，为异源酶设置一个简化的蛋白成本：

```text
单步反应的酶成本系数 = MW / (kcat × 3600)
反应通量 × 酶成本系数 = 该反应占用的近似蛋白量
```

其中 `MW` 使用 `g/mmol`（数值上等于 `kDa`），`kcat` 使用 `s⁻¹`，乘以 `3600` 是把秒换算成小时。按当前实现，蛋白池数值可以近似理解为 `g enzyme/gDW`，但它仍然是简化伪池，不等同于实测蛋白组数据。

| 参数 | 默认值 | 用来做什么 | 默认值来源 |
|---|---|---|---|
| `--heterologous-pools` | `0.005,0.01,0.02,0.05` | 扫描全部异源酶可使用的总蛋白预算 | 第一轮敏感性分析范围；用于观察蛋白预算是否成为瓶颈，没有对应的本项目蛋白组实测值 |
| `--p450-pools` | `0.001,0.002,0.005,0.01` | 单独限制 P450 类膜蛋白可使用的预算 | 第一轮敏感性分析范围；因为 P450 表达、膜定位和电子传递通常更难，所以另外设置一个更小的池，但具体数值仍需实验校准 |
| `--heterologous-pool` | `0.02` | 模块扫描时使用的固定总异源蛋白池 | 从上述扫描范围中选取的中间代表值，不是实测蛋白含量 |
| `--p450-pool` | `0.005` | 模块扫描时使用的固定 P450 池 | 从上述扫描范围中选取的中间代表值，不是实测膜蛋白含量 |
| `--expression-atp-cost-per-rel` | `0.02 mmol/gDW/h`/相对表达单位 | 把翻译、折叠和维持压力近似为固定 ATP 水解通量 | 项目中的试探性负担系数。当前没有表达量和 ATP 成本数据，因此只能用于有/无负担或不同系数之间的比较 |
| `--naringenin-feed` | `1.0 mmol/gDW/h` | 给后段 FNS/F3'H 模块提供柚皮素，使它能脱离上游模块单独测试 | 诊断性供给上限，目的是避免底物完全缺失；不是培养基中柚皮素浓度，也不是实测摄取速率 |
| `endpoint_min` | `0.001 mmol/gDW/h` | 要求某个模块至少生成少量中间产物或木犀草素 | `data/module_design.csv` 中设置的低水平测试通量，用来确认模块是否连通和可行 |
| `--enzymes` | `data/enzymes.csv` | 指定酶参数表 | 项目输入文件，可在获得可靠酶学数据后替换 |
| `--modules` | `data/module_design.csv` | 指定每个模块启用哪些反应 | 按本项目的实验拆分思路设置，包括 WT、底盘、前体、苯丙烷、前段黄酮、后段黄酮和完整通路 |

### 6.3 `enzymes.csv` 中各列的含义

| 列名 | 作用 |
|---|---|
| `reaction_id` | 与模型中的异源反应编号对应 |
| `gene` | 便于阅读的基因或酶名称 |
| `mw_kda` | 有效分子量，用于把酶用量换算成蛋白质量成本 |
| `kcat_s` | 每秒周转数；数值越低，同样通量需要的酶越多 |
| `is_p450` | 为 `1` 时，同时消耗总异源蛋白池和 P450 蛋白池 |
| `module` | 该酶所属的通路模块，仅用于整理和解释 |
| `relative_expression` | 相对表达权重，只用于计算固定 ATP 表达成本，不是转录组或蛋白组测量值 |

当前表中的数值如下：

| 酶/反应 | MW (kDa) | `kcat` (s⁻¹) | P450 | 相对表达权重 |
|---|---:|---:|---:|---:|
| `PheA_fbr` | 43 | 20 | 0 | 0.80 |
| `PAL_Li` | 77 | 10 | 0 | 0.60 |
| `C4H_ATR2` | 120 | 0.5 | 1 | 0.50 |
| `4CL_Morus` | 60 | 5 | 0 | 0.70 |
| `CHS_Carya` | 43 | 5 | 0 | 1.00 |
| `CHI_virus` | 25 | 20 | 0 | 0.80 |
| `FNS_Sm` | 58 | 0.5 | 1 | 1.00 |
| `F3pH_ATR2_main` | 120 | 0.3 | 1 | 0.40 |

这些 `MW`、`kcat` 和相对表达量是第一版建模时使用的有效值或保守估计。目前仓库中没有保存它们逐项对应的序列计算记录、BRENDA 条目或原始实验文献，因此不能写成“已经由文献证实的酶学参数”。其中 `C4H_ATR2` 和 `F3pH_ATR2_main` 的 `120 kDa` 还是把 P450 与伴随还原酶负担合并考虑的有效值。后续应优先用实际构建所用的蛋白序列重新计算 MW，并用同物种、同底物且条件接近的 `kcat` 或实测酶活替换。

### 6.4 代码内部的情景参数

下面这些参数没有暴露成常用命令行选项，但会影响特定比较：

| 参数 | 当前设置 | 作用和来源 |
|---|---:|---|
| 异源反应上限 | PheA 8、PAL 6、C4H 5、4CL 7、CHS 10、CHI 8、FNS 10、F3'H 4 | 早期项目中的占位容量上限。低通量扫描时一般不会首先达到这些上限，但做最大产能分析时可能影响结果；目前没有逐项实验来源 |
| 转运/交换上限 | `1000` | 常见的“近似不限制”大数值，只表示不希望该反应成为瓶颈，不代表真实转运能力 |
| F3'H 旁路上限 | `0.8` | 将柚皮素到圣草酚的旁路设得低于主反应，用于情景比较；属于项目假设 |
| 容量调整倍数 | ARO3/4：1.5；TKL1：2.0；TAL1：2.0；PYK1/2：0.6；PCK1：1.5；ACC1：1.5 | 来自早期改造方案的方向性设定，用于比较增强或削弱后的趋势，不是蛋白表达倍数或实测酶活倍数 |
| 酪氨酸摄取救援 | `1.0 mmol/gDW/h` | 判断旧 PheA 一步式模型中的不生长是否可被外源酪氨酸救回，属于诊断性对照 |
| ARO7 弱化场景 | WT 模型通量的 10% | 用来比较完全敲除与部分保留活性的差别，是敏感性测试值 |
| TYR1 25% 上限场景 | WT TYR1 通量的 25% | 25% 来自项目报告对 PheA 文献 Fig. 4 中“HMA 约降至四分之一”的读图总结，只用于测试“若把它错误等同为 TYR1 上限会怎样”；这不是文献直接给出的 TYR1 参数，HMA 也不能直接当作维持生长的 TYR1 通量 |
| PheA CM/PDT 上限 | `1000` | 把双功能 PheA 的两步反应设为近似不受上限限制，用于先验证计量可行性；不是 PheA 的真实催化能力 |

### 6.5 修改参数的示例

大部分运行脚本都可以用命令行参数修改条件。例如：

```powershell
python src/run_fba.py --luteolin-demands 0,0.001,0.01
python src/run_ecfba_scan.py --heterologous-pools 0.01,0.02 --p450-pools 0.002,0.005
python src/run_module_scan.py --heterologous-pool 0.02 --p450-pool 0.005
```

可以用 `-h` 查看某个脚本支持的参数：

```powershell
python src/run_ecfba_scan.py -h
```

修改时建议一次只改变一类参数，并保留 WT 和原默认条件作为对照。若后续获得实验数据，建议按“真实葡萄糖/氧气摄取速率 → 实际蛋白 MW → 酶学 `kcat` → 蛋白组或表达量 → 生长和产量”的顺序逐步替换假设值。

## 7. 默认参数和建模方法的参考来源

1. yeastGEM v9.0.2 模型版本：[SysBioChalmers/yeast-GEM v9.0.2](https://github.com/SysBioChalmers/yeast-GEM/releases/tag/v9.0.2)。本项目实际使用的版本也可以在 `models/yeastGEM.xml` 第 3 行核对。
2. yeast-GEM / Yeast9 模型说明：Zhang C. 等，*Yeast9: a consensus genome-scale metabolic model for S. cerevisiae curated by the community*，2024，[DOI: 10.1038/s44320-024-00060-7](https://doi.org/10.1038/s44320-024-00060-7)。
3. FBA 求解和反应上下界的实现：[COBRApy 官方文档](https://opencobra.github.io/cobrapy/)。
4. 酶容量约束的基本思想：Sánchez B. J. 等，*Improving the phenotype predictions of a yeast genome-scale metabolic model by incorporating enzymatic constraints*，2017，[PubMed 28779005](https://pubmed.ncbi.nlm.nih.gov/28779005/)。
5. 完整 GECKO 的建模和校准流程：Chen Y. 等，*Reconstruction, simulation and analysis of enzyme-constrained metabolic models using GECKO Toolbox 3.0*，2024，[DOI: 10.1038/s41596-023-00931-7](https://doi.org/10.1038/s41596-023-00931-7)。
6. `PheA^fbr` 双功能、`aro7Δ` 生长及 MA/HMA 对照依据：Reifenrath M. 等，*Bacterial bifunctional chorismate mutase-prephenate dehydratase PheA increases flux into the yeast phenylalanine pathway and improves mandelic acid production*，2018，[DOI: 10.1016/j.mec.2018.e00079](https://doi.org/10.1016/j.mec.2018.e00079)。
