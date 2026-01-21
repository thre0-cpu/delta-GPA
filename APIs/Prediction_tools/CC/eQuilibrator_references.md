# 📚 eQuilibrator 参考文档

根据 eQuilibrator 官方文档和 Python API 文档最佳实践，由Claude 4.5 Sonnet编写

---

## 目录

1. [快速开始](#快速开始)
2. [核心模块](#核心模块)
3. [模型类](#模型类)
4. [化合物与反应](#化合物与反应)
5. [工具函数](#工具函数)
6. [常量与默认值](#常量与默认值)
7. [完整使用示例](#完整使用示例)
---

## 快速开始

eQuilibrator API 是一个用于计算生化反应热力学参数的 Python 包，支持标准吉布斯自由能、平衡常数等计算 [3]。

### 基本安装

```bash
pip install equilibrator-api
```

### 简单示例

```python
from equilibrator_api import ComponentContribution

# 初始化
cc = ComponentContribution()

# 搜索化合物
glucose = cc.search_compound("glucose")

# 解析反应
rxn = cc.parse_reaction_formula("glucose + ATP => glucose-6-phosphate + ADP")

# 计算标准吉布斯自由能
dg_prime = cc.standard_dg_prime(rxn)
print(f"ΔG'° = {dg_prime}")
```

---

## 核心模块

### 1. `equilibrator_api.component_contribution`

Component Contribution 方法的包装器，用于预测吉布斯能量 [4]。

#### 主类：`ComponentContribution`

**继承自:** `object`

**描述:** GibbsEnergyPredictor 的包装类，包含不同相态化合物的默认条件。

##### 构造函数

```python
ComponentContribution(
    rmse_inf: Q_ = default_rmse_inf,
    ccache: Optional[CompoundCache] = None,
    predictor: Optional[GibbsEnergyPredictor] = None
)
```

**参数:**
- `rmse_inf` (Quantity): 误差协方差矩阵的乘数因子，用于 Component Contribution 范围外的反应（默认: 1e-5 kJ/mol）
- `ccache` (CompoundCache, optional): 化合物缓存对象
- `predictor` (GibbsEnergyPredictor, optional): 自定义预测器

##### 属性（Properties）

| 属性 | 类型 | 描述 |
|------|------|------|
| `p_h` | `Q_` | 获取 pH 值 |
| `p_mg` | `Q_` | 获取 pMg 值 |
| `ionic_strength` | `Q_` | 获取离子强度 |
| `temperature` | `Q_` | 获取温度 |
| `RT` | `Q_` | 获取 RT 值（气体常数×温度）|

##### 核心方法

###### 化合物搜索

```python
def get_compound(compound_id: str) -> Union[Compound, None]:
    """
    使用数据库命名空间和登录号获取化合物
    
    参数:
        compound_id: 化合物标识符（如 "KEGG:C00002"）
    
    返回:
        Compound 对象，如果未找到则返回 None
    
    示例:
        >>> cc = ComponentContribution()
        >>> atp = cc.get_compound("KEGG:C00002")
        >>> print(atp.formula)
    """
```

```python
def search_compound(query: str) -> Union[None, Compound]:
    """
    尝试找到最匹配名称的化合物
    
    参数:
        query: 化合物名称（近似匹配）
    
    返回:
        最佳匹配的 Compound 对象
    
    示例:
        >>> cc = ComponentContribution()
        >>> glucose = cc.search_compound("glucose")
        >>> print(glucose.id)
    """
```

```python
def get_compound_by_inchi(inchi: str) -> Union[Compound, None]:
    """使用 InChI 获取化合物"""
```

```python
def search_compound_by_inchi_key(inchi_key: str) -> List[Compound]:
    """使用 InChI Key 搜索化合物"""
```

###### 反应解析

```python
def parse_reaction_formula(formula: str) -> PhasedReaction:
    """
    使用精确匹配解析反应文本
    
    参数:
        formula: 包含反应式的字符串
    
    返回:
        PhasedReaction 对象
    
    示例:
        >>> cc = ComponentContribution()
        >>> rxn = cc.parse_reaction_formula("ATP + H2O => ADP + Pi")
    """
```

```python
def search_reaction(formula: str) -> PhasedReaction:
    """
    使用化合物名称搜索反应（近似匹配）
    
    参数:
        formula: 包含反应式的字符串
    
    返回:
        PhasedReaction 对象
    """
```

###### 热力学计算

```python
def standard_dg(reaction: PhasedReaction) -> ureg.Measurement:
    """
    计算反应的化学反应能
    
    返回:
        标准吉布斯自由能 ΔG° (kJ/mol) 及标准误差
        95% 置信区间 = ±1.96 × 标准误差
    
    示例:
        >>> cc = ComponentContribution()
        >>> rxn = cc.parse_reaction_formula("glucose => 2 pyruvate")
        >>> dg0 = cc.standard_dg(rxn)
        >>> print(f"ΔG° = {dg0.value:.2f} ± {dg0.error:.2f} kJ/mol")
    """
```

```python
def standard_dg_prime(reaction: PhasedReaction) -> ureg.Measurement:
    """
    计算反应的转化反应能
    
    返回:
        标准转化吉布斯自由能 ΔG'° (kJ/mol) 及标准误差
        考虑了 pH、离子强度等因素
    """
```

```python
def physiological_dg_prime(reaction: PhasedReaction) -> ureg.Measurement:
    """
    计算单个反应的生理 ΔG'm
    
    假设:
        - 所有水相反应物浓度为 1 mM
        - 气相反应物压力为 1 mbar
        - 其他反应物为标准浓度
    
    返回:
        生理条件下的 ΔG' 及标准误差
    """
```

```python
def dg_prime(reaction: PhasedReaction) -> ureg.Measurement:
    """
    计算单个反应的 ΔG'
    
    考虑实际浓度/压力条件
    """
```

```python
def standard_dg_multi(
    reactions: List[PhasedReaction],
    uncertainty_representation: str = 'cov'
) -> Tuple[np.ndarray, np.ndarray]:
    """
    计算多个反应的化学反应能
    
    参数:
        reactions: PhasedReaction 对象列表
        uncertainty_representation: 不确定性表示方式
            - 'cov': 完整协方差矩阵
            - 'precision': 精度矩阵（协方差矩阵的逆）
            - 'sqrt': 协方差的平方根
            - 'fullrank': 满秩平方根（压缩形式）
    
    返回:
        (standard_dg, dg_uncertainty): 标准反应能数组和不确定性矩阵
    """
```

```python
def standard_dg_prime_multi(
    reactions: List[PhasedReaction],
    uncertainty_representation: str = 'cov',
    minimize_norm: bool = False
) -> Tuple[Q_, Q_]:
    """
    计算多个反应的转化反应能
    
    参数:
        reactions: PhasedReaction 对象列表
        uncertainty_representation: 不确定性表示方式
        minimize_norm: 是否使用正交投影最小化结果向量的范数
    
    返回:
        (standard_dg_prime, dg_uncertainty)
    """
```

###### 电化学计算

```python
def standard_e_prime(reaction: PhasedReaction) -> ureg.Measurement:
    """
    计算单个半反应的标准电势 E'°
    
    返回:
        标准电势 (V) 及标准误差
        95% 置信区间 = ±1.96 × 标准误差
    """
```

```python
def physiological_e_prime(reaction: PhasedReaction) -> ureg.Measurement:
    """计算生理条件下的电势 E'm"""
```

```python
def e_prime(reaction: PhasedReaction) -> ureg.Measurement:
    """计算实际条件下的电势 E'"""
```

###### 敏感性分析

```python
def dgf_prime_sensitivity_to_p_h(compound: Compound) -> ureg.Quantity:
    """
    计算化学生成能对 pH 的敏感性
    
    返回:
        ∂(ΔGf')/∂pH 的导数 (kJ/mol)
    """
```

```python
def dg_prime_sensitivity_to_p_h(reaction: PhasedReaction) -> ureg.Quantity:
    """
    计算化学反应能对 pH 的敏感性
    
    返回:
        ∂(ΔGr')/∂pH 的导数 (kJ/mol)
    """
```

###### 可逆性指数

```python
def ln_reversibility_index(reaction: PhasedReaction) -> ureg.Measurement:
    """
    计算单个反应的可逆性指数（ln Γ）
    
    返回:
        自然对数尺度的可逆性指数
    """
```

###### 高级功能

```python
def balance_by_oxidation(reaction: PhasedReaction) -> PhasedReaction:
    """
    将不平衡反应转换为氧化反应
    
    通过向两侧添加 H2O、O2、Pi、CO2 和 NH4+ 来平衡
    """
```

```python
def get_oxidation_reaction(compound: Compound) -> PhasedReaction:
    """
    为单个化合物生成氧化反应
    
    使用 O2 生成氧化反应。对于含 N 原子的化合物，
    产物为 NH3（而非 N2）以代表生物过程
    """
```

```python
def multicompartmental_standard_dg_prime(
    reaction_inner: PhasedReaction,
    reaction_outer: PhasedReaction,
    e_potential_difference: Q_,
    p_h_outer: Q_,
    ionic_strength_outer: Q_,
    p_mg_outer: Q_ = default_physiological_p_mg,
    tolerance: float = 0.0
) -> ureg.Measurement:
    """
    计算多室反应的转化能
    
    基于 Haraldsdóttir et al. 2012 的方程
    (https://doi.org/10.1016/j.bpj.2012.02.032)
    
    参数:
        reaction_inner: 内室半反应
        reaction_outer: 外室半反应
        e_potential_difference: 外室与内室的电势差
        p_h_outer: 外室 pH
        ionic_strength_outer: 外室离子强度
        p_mg_outer: 外室 pMg（可选）
        tolerance: 识别内外反应不平衡的容差（默认 = 0）
    
    返回:
        转运反应的吉布斯自由能变化
    """
```

###### 分析工具

```python
def dg_analysis(reaction: PhasedReaction) -> List[Dict[str, object]]:
    """
    获取 Component Contribution 估算过程的分析
    
    返回:
        分析结果的字典列表
    """
```

```python
def is_using_group_contribution(reaction: PhasedReaction) -> bool:
    """
    检查是否需要基团贡献法来获取该反应的 ΔG
    
    返回:
        如果需要基团贡献法则返回 True
    """
```

```python
def standard_dg_formation(compound: Compound) -> Tuple[Optional[float], Optional[np.ndarray]]:
    """
    获取化合物生成能的 (μ, σ) 预测
    
    参数:
        compound: Compound 对象
    
    返回:
        - mu: 标准生成吉布斯能估计的均值
        - sigma_fin: 表示协方差矩阵平方根的向量（不确定性）
        - sigma_inf: 表示协方差矩阵无限不确定性特征值的向量
    """
```

##### 静态方法

```python
@staticmethod
def legacy() -> ComponentContribution:
    """
    使用旧版本初始化 ComponentContribution 对象
    
    旧版本用于与 equilibrator-api 旧版本（0.2.x - 0.3.1）兼容
    从 0.3.2 开始，由于改进的 Mg2+ 浓度模型，预测有显著变化
    
    返回:
        ComponentContribution 对象
    """
```

```python
@staticmethod
def initialize_custom_version(
    rmse_inf: Q_ = default_rmse_inf,
    ccache_settings: ZenodoSettings = DEFAULT_COMPOUND_CACHE_SETTINGS,
    cc_params_settings: ZenodoSettings = DEFAULT_CC_PARAMS_SETTINGS
) -> ComponentContribution:
    """
    使用自定义 Zenodo 版本初始化 ComponentContribution
    
    参数:
        rmse_inf: 误差协方差矩阵的乘数因子
        ccache_settings: 化合物缓存的 Zenodo 设置
        cc_params_settings: CC 参数的 Zenodo 设置
    
    返回:
        ComponentContribution 对象
    """
```

```python
@staticmethod
def parse_formula_side(s: str) -> Dict[str, float]:
    """解析反应式的一侧"""
```

```python
@staticmethod
def parse_formula(formula: str) -> Dict[str, float]:
    """解析双侧反应式"""
```

```python
@staticmethod
def create_stoichiometric_matrix_from_reaction_formulas(
    formulas: Iterable[str]
) -> pd.DataFrame:
    """
    构建化学计量矩阵
    
    参数:
        formulas: 反应式的字符串表示
    
    返回:
        DataFrame，索引为化合物 ID，列为反应 ID
    """
```

```python
def create_stoichiometric_matrix_from_reaction_objects(
    reactions: Iterable[PhasedReaction]
) -> pd.DataFrame:
    """
    从反应对象构建化学计量矩阵
    
    参数:
        reactions: PhasedReaction 对象集合
    
    返回:
        DataFrame，索引为化合物，列为反应
    """
```

##### 辅助函数

```python
def find_most_abundant_ms(
    cpd: Compound,
    p_h: Q_,
    p_mg: Q_,
    ionic_strength: Q_,
    temperature: Q_
) -> CompoundMicrospecies:
    """基于转化能找到最丰富的微物种"""
```

```python
def predict_protons_and_charge(
    rxn: PhasedReaction,
    p_h: Q_,
    p_mg: Q_,
    ionic_strength: Q_,
    temperature: Q_
) -> Tuple[float, float, float]:
    """找到转运半反应的质子数和电荷"""
```

---

### 2. `equilibrator_api.model.bounds`

定义化合物的上下界限。

#### 类：`BaseBounds`

**继承自:** `object`

**描述:** 声明边界的基类。

##### 抽象方法

```python
@abstractmethod
def copy(self) -> BaseBounds:
    """返回 self 的深拷贝"""
```

```python
@abstractmethod
def get_lower_bound(compound: Union[str, Compound]) -> Q_:
    """
    获取该键的下界
    
    参数:
        compound: 化合物对象或字符串
    """
```

```python
@abstractmethod
def get_upper_bound(compound: Union[str, Compound]) -> Q_:
    """
    获取该键的上界
    
    参数:
        compound: 化合物对象或字符串
    """
```

##### 实例方法

```python
def get_lower_bounds(compounds: Iterable[Union[str, Compound]]) -> Iterable[Q_]:
    """
    按顺序获取一组键的下界
    
    参数:
        compounds: Compound 或字符串的可迭代对象
    
    返回:
        下界的可迭代对象
    """
```

```python
def get_upper_bounds(compounds: Iterable[Union[str, Compound]]) -> Iterable[Q_]:
    """按顺序获取一组键的上界"""
```

```python
def get_bound_tuple(compound: Union[str, Compound]) -> Tuple[Q_, Q_]:
    """
    获取该键的上下界
    
    参数:
        compound: Compound 对象或字符串
    
    返回:
        二元组 (下界, 上界)
    """
```

```python
def get_bounds(compounds: Iterable[Union[str, Compound]]) -> Tuple[Iterable[Q_], Iterable[Q_]]:
    """
    获取一组化合物的边界
    
    参数:
        compounds: Compound 的可迭代对象
    
    返回:
        二元组 (下界列表, 上界列表)
    """
```

```python
@staticmethod
def conc2ln_conc(b: Q_) -> float:
    """
    将浓度转换为对数浓度
    
    参数:
        b: 浓度
    
    返回:
        对数浓度
    """
```

```python
def get_ln_bounds(compounds: Iterable[Union[str, Compound]]) -> Tuple[Iterable[float], Iterable[float]]:
    """
    获取一组化合物的对数边界
    
    返回:
        二元组 (对数下界, 对数上界)
    """
```

```python
def get_ln_lower_bounds(compounds: Iterable[Union[str, Compound]]) -> Iterable[float]:
    """获取对数下界"""
```

```python
def get_ln_upper_bounds(compounds: Iterable[Union[str, Compound]]) -> Iterable[float]:
    """获取对数上界"""
```

```python
def set_bounds(compound: Union[str, Compound], lb: Q_, ub: Q_) -> None:
    """
    设置特定键的边界
    
    参数:
        compound: Compound 或字符串
        lb: 下界值
        ub: 上界值
    """
```

#### 类：`Bounds`

**继承自:** `BaseBounds`

**描述:** 包含各种键的上下界，允许默认值。

##### 构造函数

```python
Bounds(
    lower_bounds: Dict[Union[str, Compound], Q_] = None,
    upper_bounds: Dict[Union[str, Compound], Q_] = None,
    default_lb: Q_ = default_conc_lb,
    default_ub: Q_ = default_conc_ub
)
```

##### 类属性

```python
DEFAULT_BOUNDS  # 默认边界
```

##### 类方法

```python
@classmethod
def from_csv(
    f: TextIO,
    comp_contrib: ComponentContribution,
    default_lb: Q_ = default_conc_lb,
    default_ub: Q_ = default_conc_ub
) -> Bounds:
    """
    从 CSV 文件读取 Bounds
    
    参数:
        f: 打开的 .csv 文件流
        comp_contrib: 用于解析化合物登录号
        default_lb: 默认下界
        default_ub: 默认上界
    """
```

##### 实例方法

```python
def to_data_frame() -> pd.DataFrame:
    """将边界列表转换为 Pandas DataFrame"""
```

```python
def check_bounds() -> None:
    """断言边界有效（即 lb <= ub）"""
```

```python
def copy() -> Bounds:
    """返回 self 的深拷贝"""
```

```python
def get_lower_bound(compound: Union[str, Compound]) -> Q_:
    """获取该化合物的下界"""
```

```python
def get_upper_bound(compound: Union[str, Compound]) -> Q_:
    """获取该化合物的上界"""
```

##### 静态方法

```python
@staticmethod
def get_default_bounds(comp_contrib: ComponentContribution) -> Bounds:
    """
    返回预定列表的默认上下界
    
    参数:
        comp_contrib: ComponentContribution 对象
    
    返回:
        具有默认值的 Bounds 对象
    """
```

---

### 3. `equilibrator_api.model.model`

带热力学的基本化学计量模型。

#### 类：`StoichiometricModel`

**继承自:** `object`

**描述:** 带热力学的基本化学计量模型，设计为 'Pathway' 的基础模型，后者还包括通量方向和大小。

##### 构造函数

```python
StoichiometricModel(
    S: pd.DataFrame,
    compound_dict: Dict[str, Compound],
    reaction_dict: Dict[str, Reaction],
    comp_contrib: Optional[ComponentContribution] = None,
    standard_dg_primes: Optional[Q_] = None,
    dg_sigma: Optional[Q_] = None,
    bounds: Optional[Bounds] = None,
    config_dict: Optional[Dict[str, str]] = None
)
```

##### 类属性

```python
MINIMAL_STDEV = 0.001
```

##### 实例方法

```python
def configure() -> None:
    """配置 Component Contribution 水相条件"""
```

##### 属性（Properties）

```python
@property
def compound_ids() -> Iterable[str]:
    """获取化合物 ID 列表"""
```

```python
@property
def compounds() -> Iterable[Compound]:
    """获取 Compound 对象列表"""
```

```python
@property
def compound_df() -> pd.DataFrame:
    """
    获取包含所有化合物数据的 DataFrame
    
    列:
        - compound_id
        - lower_bound
        - upper_bound
    """
```

```python
@property
def reaction_ids() -> Iterable[str]:
    """获取反应 ID 列表"""
```

```python
@property
def reactions() -> Iterable[Reaction]:
    """获取 Reaction 对象列表"""
```

```python
@property
def reaction_formulas() -> Iterable[str]:
    """
    迭代所有反应式
    
    返回:
        反应式字符串
    """
```

```python
@property
def reaction_df() -> pd.DataFrame:
    """
    获取包含所有反应数据的 DataFrame
    
    列:
        - reaction_id
        - reaction_formula
        - standard_dg_prime
    """
```

```python
@property
def bounds() -> Tuple[Iterable[Q_], Iterable[Q_]]:
    """
    获取浓度边界
    
    化合物顺序按化学计量矩阵索引
    
    返回:
        (下界, 上界) 元组
    """
```

```python
@property
def bound_df() -> pd.DataFrame:
    """获取包含所有边界数据的 DataFrame"""
```

```python
@property
def ln_conc_lb() -> np.array:
    """
    获取浓度的对数下界
    
    化合物顺序按化学计量矩阵索引
    
    返回:
        对数下界的 NumPy 数组
    """
```

```python
@property
def ln_conc_ub() -> np.ndarray:
    """
    获取浓度的对数上界
    
    返回:
        对数上界的 NumPy 数组
    """
```

```python
@property
def ln_conc_mu() -> np.array:
    """
    基于边界获取对数浓度分布的均值
    
    返回:
        对数浓度均值的 NumPy 数组
    """
```

```python
@property
def ln_conc_sigma() -> np.array:
    """
    基于边界获取对数浓度分布的标准差
    
    返回:
        对数浓度标准差的 NumPy 数组
    """
```

##### 实例方法

```python
def update_standard_dgs() -> None:
    """
    计算标准 G' 值和不确定性
    
    使用 Component Contribution 方法
    """
```

```python
def set_bounds(cid: str, lb: Optional[Q_] = None, ub: Optional[Q_] = None) -> None:
    """
    设置化合物的上下界
    
    参数:
        cid: 化合物 ID
        lb: 新的浓度下界（如果为 None 则忽略）
        ub: 新的浓度上界（如果为 None 则忽略）
    """
```

```python
def get_bounds(cid: str) -> Tuple[Q_, Q_]:
    """
    获取化合物的上下界
    
    参数:
        cid: 化合物 ID
    
    返回:
        (下界, 上界)
    """
```

##### 静态方法

```python
@staticmethod
def read_thermodynamics(thermo_sbtab: SBtabTable, config_dict: Dict[str, str]) -> Dict[str, Q_]:
    """
    从 SBtab 读取 'thermodynamics' 表
    
    参数:
        thermo_sbtab: 包含热力学数据的 SBtabTable
        config_dict: 包含配置参数的字典
    
    返回:
        将反应 ID 映射到标准 G' 值的字典
    """
```

##### 类方法

```python
@classmethod
def from_network_sbtab(
    filename: Union[str, SBtabDocument],
    comp_contrib: Optional[ComponentContribution] = None,
    freetext: bool = True,
    bounds: Optional[Bounds] = None
) -> StoichiometricModel:
    """
    使用仅包含 'network' 的 SBtab 初始化 Pathway 对象
    
    参数:
        filename: 包含 SBtabDocument 的文件名（或 SBtabDocument 对象本身）
        comp_contrib: 用于解析和搜索反应的 ComponentContribution 对象
        freetext: 标志，指示反应是否以自由文本给出（即化合物的通用名称）
                  或通过标准数据库登录号（默认: True）
        bounds: 代谢物浓度边界（默认使用 equilibrator-api 中的 "data/cofactors.csv" 文件）
    
    返回:
        Pathway 对象
    """
```

```python
@classmethod
def from_sbtab(
    filename: Union[str, SBtabDocument],
    comp_contrib: Optional[ComponentContribution] = None
) -> StoichiometricModel:
    """
    解析 SBtabDocument 并返回 StoichiometricModel
    
    参数:
        filename: 包含 SBtabDocument 的文件名（或 SBtabDocument 对象本身）
        comp_contrib: 用于解析和搜索反应的 ComponentContribution 对象
    
    返回:
        基于配置 SBtab 的 StoichiometricModel 对象
    """
```

##### 导出方法

```python
def to_sbtab() -> SBtabDocument:
    """将模型导出为 SBtabDocument"""
```

```python
def write_sbtab(filename: str) -> None:
    """将路径写入 SBtab 文件"""
```

---

## 化合物与反应

### 4. `equilibrator_api.phased_compound`

从 equilibrator_cache.models.compound.Compound 继承并添加相态。

#### 常量

```python
AQUEOUS_PHASE_NAME = "aqueous"
GAS_PHASE_NAME = "gas"
LIQUID_PHASE_NAME = "liquid"
SOLID_PHASE_NAME = "solid"
REDOX_PHASE_NAME = "redox"
```

#### 类型别名

```python
PhaseInfo
MicroSpecie
PHASE_INFO_DICT
NON_AQUEOUS_COMPOUND_DICT
PHASED_COMPOUND_DICT
CARBONATE_INCHIS
```

#### 类：`Condition`

**继承自:** `object`

**描述:** 定义化合物条件的类，即相态和丰度。

##### 构造函数

```python
Condition(phase: str, abundance: ureg.Quantity = None)
```

##### 属性

```python
@property
def phase() -> str:
    """返回相态"""
```

```python
@property
def abundance() -> ureg.Quantity:
    """返回丰度"""
```

```python
@property
def standard_abundance() -> ureg.Quantity:
    """返回该相态的标准丰度"""
```

```python
@property
def physiological_abundance() -> ureg.Quantity:
    """返回该相态的默认生理丰度"""
```

```python
@property
def dimensionality() -> str:
    """
    返回该相态丰度的量纲
    
    例如水相为 [concentration]，气相为 [pressure]
    
    返回:
        该相态的量纲，如果丰度固定则返回 None
    """
```

```python
@property
def ln_abundance() -> float:
    """返回给定丰度与标准丰度比值的对数"""
```

```python
@property
def ln_physiological_abundance() -> float:
    """返回生理丰度与标准丰度比值的对数"""
```

```python
@property
def is_physiological() -> bool:
    """
    如果丰度与生理丰度相同则返回 True
    
    返回:
        如果丰度处于生理条件，或该相态丰度固定，则返回 True
    """
```

##### 方法

```python
def reset_abundance() -> None:
    """将丰度重置为标准丰度"""
```

#### 类：`PhasedCompound`

**继承自:** `object`

**描述:** 结合 equilibrator_api Compound 和 Condition 的类。

##### 构造函数

```python
PhasedCompound(compound: Compound, condition: Condition = None)
```

##### 静态方法

```python
@staticmethod
def get_default(compound: Compound) -> Condition:
    """
    获取化合物的默认相态
    
    参数:
        compound: Compound 对象
    
    返回:
        默认相态
    """
```

##### 属性

```python
@property
def atom_bag() -> dict:
    """获取化合物的原子袋"""
```

```python
@property
def smiles() -> str:
    """获取化合物的 SMILES"""
```

```python
@property
def inchi() -> str:
    """获取化合物的 InChI"""
```

```python
@property
def inchi_key() -> str:
    """获取化合物的 InChIKey"""
```

```python
@property
def id() -> int:
    """获取化合物的 equilibrator 内部 ID"""
```

```python
@property
def formula() -> str:
    """获取化学式"""
```

```python
@property
def mass() -> float:
    """获取化学分子质量"""
```

```python
@property
def phase() -> str:
    """获取相态"""
```

```python
@property
def html_formula() -> str:
    """获取 HTML 格式的化学式"""
```

```python
@property
def phase_shorthand() -> str:
    """获取相态简写（例如液相为 'l'）"""
```

```python
@property
def possible_phases() -> Tuple[str]:
    """获取该化合物的可能相态"""
```

```python
@property
def abundance() -> ureg.Quantity:
    """获取丰度"""
```

```python
@property
def ln_abundance() -> float:
    """返回丰度的对数（用于热力学计算）"""
```

```python
@property
def ln_physiological_abundance() -> float:
    """返回默认生理丰度的对数"""
```

```python
@property
def is_physiological() -> bool:
    """检查丰度是否为生理条件"""
```

##### 方法

```python
def get_stored_standard_dgf_prime(
    p_h: ureg.Quantity,
    ionic_strength: ureg.Quantity,
    temperature: ureg.Quantity,
    p_mg: ureg.Quantity
) -> ureg.Quantity:
    """
    返回该相态化合物的存储生成能
    
    仅在存在时返回，否则返回 None
    （稍后将使用 component-contribution 获取反应能）
    
    参数:
        p_h: pH
        ionic_strength: 离子强度
        temperature: 温度
        p_mg: pMg
    
    返回:
        standard_dgf_prime (kJ/mol)
    """
```

```python
def get_stored_standard_dgf() -> ureg.Quantity:
    """
    返回该相态化合物的存储生成能
    
    返回:
        standard_dgf (kJ/mol)
    """
```

```python
def get_stored_microspecie() -> MicroSpecie:
    """
    获取存储的微物种（来自 PHASED_COMPOUND_DICT）
    
    返回:
        包含存储生成能的 MicroSpecie 命名元组，
        如果该化合物在此相态没有存储值则返回 None
    """
```

```python
def serialize() -> dict:
    """
    返回所有化合物热力学数据的序列化版本
    
    返回:
        包含所有微物种数据的字典列表
    """
```

#### 类：`Proton`

**继承自:** `PhasedCompound`

**描述:** 专门用于质子的类。

##### 属性

```python
@property
def abundance() -> ureg.Quantity:
    """获取丰度"""
```

```python
@property
def ln_physiological_abundance() -> float:
    """返回默认生理丰度的对数"""
```

```python
@property
def ln_abundance() -> float:
    """返回丰度的对数（用于热力学计算）"""
```

#### 类：`RedoxCarrier`

**继承自:** `PhasedCompound`

**描述:** 专门用于氧化还原载体（具有给定电势）的类。

##### 构造函数

```python
RedoxCarrier(compound: Compound, potential: Optional[ureg.Quantity] = None)
```

##### 方法

```python
def get_stored_standard_dgf_prime(...) -> ureg.Quantity:
    """获取标准生成 G'"""
```

```python
def get_stored_standard_dgf() -> ureg.Quantity:
    """获取标准生成 G"""
```

##### 属性

```python
@property
def atom_bag() -> dict:
    """获取化合物的原子袋"""
```

```python
@property
def ln_abundance() -> float:
    """返回丰度的对数（用于热力学计算）"""
```

```python
@property
def ln_physiological_abundance() -> float:
    """返回默认生理丰度的对数"""
```

```python
@property
def is_physiological() -> bool:
    """检查丰度是否为生理条件"""
```

---

### 5. `equilibrator_api.phased_reaction`

从 equilibrator_cache.reaction.Reaction 继承并添加相态。

#### 类：`PhasedReaction`

**继承自:** `equilibrator_cache.Reaction`

**描述:** Reaction 的子类，添加了相态和丰度。

##### 类属性

```python
REACTION_COUNTER = 0
```

##### 构造函数

```python
PhasedReaction(
    sparse: Dict[Compound, float],
    arrow: str = '<=>',
    rid: str = None,
    sparse_with_phases: Dict[PhasedCompound, float] = None
)
```

##### 静态方法

```python
@staticmethod
def to_phased_compound(cpd: Compound) -> PhasedCompound:
    """将 Compound 对象转换为 PhasedCompound"""
```

##### 实例方法

```python
def clone() -> PhasedReaction:
    """克隆此反应对象"""
```

```python
def reverse() -> PhasedReaction:
    """返回逆反应的 PhasedReaction"""
```

```python
def get_element_data_frame() -> pd.DataFrame:
    """
    列出所有反应物的元素组成
    
    返回:
        DataFrame，列为化合物，索引为原子元素
    """
```

```python
def hash_md5(reversible: bool = True) -> str:
    """
    返回 PhasedReaction 的 MD5 哈希
    
    此哈希用于查找具有完全相同化学计量的反应
    基于 Compound ID 和系数创建唯一的公式字符串
    
    参数:
        reversible: 标志，指示反应方向是否重要
                    如果为 True，正向和反向版本将返回相同值
    
    返回:
        表示 Reaction 的唯一哈希字符串
    """
```

```python
def set_abundance(compound: Compound, abundance: ureg.Quantity):
    """设置化合物的丰度"""
```

```python
def reset_abundances():
    """将丰度重置为标准水平"""
```

```python
def set_phase(compound: Compound, phase: str):
    """设置化合物的相态"""
```

```python
def get_phased_compound(compound: Compound) -> Tuple[PhasedCompound, float]:
    """通过 Compound 对象获取 PhasedCompound 对象"""
```

```python
def get_phase(compound: Compound) -> str:
    """获取化合物的相态"""
```

```python
def get_abundance(compound: Compound) -> ureg.Quantity:
    """获取化合物的丰度"""
```

```python
def get_stoichiometry(compound: Compound) -> float:
    """获取化合物的化学计量系数"""
```

```python
def add_stoichiometry(compound: Compound, coeff: float) -> None:
    """
    添加化合物的化学计量系数
    
    如果该化合物不在反应中，则添加它
    """
```

##### 属性

```python
@property
def is_physiological() -> bool:
    """
    检查所有浓度是否为生理条件
    
    eQuilibrator 使用此函数判断是否显示调整后的 ΔG'
    （因为生理 ΔG' 总是显示，重复显示会冗余）
    
    返回:
        如果所有化合物都处于生理丰度则返回 True
    """
```

##### 高级方法

```python
def separate_stored_dg_prime(
    p_h: ureg.Quantity,
    ionic_strength: ureg.Quantity,
    temperature: ureg.Quantity,
    p_mg: ureg.Quantity
) -> Tuple[Reaction, ureg.Quantity]:
    """
    将 PhasedReaction 分为水相和其他所有相
    
    参数:
        p_h: pH
        ionic_strength: 离子强度
        temperature: 温度
        p_mg: pMg
    
    返回:
        (residual_reaction, stored_dg_prime)
        - residual_reaction: Reaction 对象（不包括具有存储值的化合物）
        - stored_dg_prime: 具有存储值的化合物的总 ΔG' (kJ/mol)
    """
```

```python
def separate_stored_dg() -> Tuple[Reaction, ureg.Quantity]:
    """
    将 PhasedReaction 分为水相和其他所有相
    
    返回:
        (residual_reaction, stored_dg)
        - residual_reaction: Reaction 对象（不包括具有存储值的化合物）
        - stored_dg: 具有存储值的化合物的总 ΔG (kJ/mol)
    """
```

```python
def dg_correction() -> ureg.Quantity:
    """
    计算反应 ΔG' 的浓度调整
    
    返回:
        以 RT 为单位的 delta G 校正
    """
```

```python
def physiological_dg_correction() -> ureg.Quantity:
    """
    计算反应 ΔG' 的浓度调整
    
    假设所有反应物处于默认生理浓度（即 1 mM）
    
    返回:
        以 RT 为单位的 delta G 校正
    """
```

```python
def serialize() -> List[dict]:
    """
    返回所有反应热力学数据的序列化版本
    
    返回:
        字典列表
    """
```

---

### 6. `equilibrator_api.reaction_parser`

反应式解析器。

#### 常量

```python
POSSIBLE_REACTION_ARROWS = [
    '<=>',  '<->',  '-->',  '<--',
    '=>',   '<=',   '->',   '<-',
    '=',    '⇄',    '⇌',    '→',    '←'
]
```

#### 函数

```python
def make_reaction_parser() -> pyparsing.Forward:
    """
    构建基于 pyparsing 的化学反应递归下降解析器
    
    返回:
        pyparsing.Forward 解析器对象
    """
```

---

### 7. `equilibrator_api.compatibility`

提供与 COBRA 兼容的函数。

#### 函数

```python
def map_cobra_reactions(
    cache: CompoundCache,
    reactions: List[cobra.Reaction],
    **kwargs
) -> Dict[str, PhasedReaction]:
    """
    将 COBRA 反应转换为 eQuilibrator 相态反应
    
    参数:
        cache: equilibrator_cache.CompoundCache 对象
        reactions: 要映射的 cobra.Reaction 列表
        kwargs: 任何其他关键字参数传递给底层代谢物映射函数
    
    返回:
        从 COBRA 反应标识符到 equilibrator 相态反应的映射字典
        （仅在可以建立映射时）
    
    另见:
        equilibrator_cache.compatibility.map_cobra_metabolites
    """
```

---

## 工具函数

### 8. `equilibrator_api.model`

#### 函数

```python
def open_sbtabdoc(filename: Union[str, SBtabDocument]) -> SBtabDocument:
    """
    将文件作为 SBtabDocument 打开
    
    检查是否已经是 SBtabDocument 对象，
    否则读取 CSV 文件并返回解析的对象
    
    参数:
        filename: 文件名或 SBtabDocument 对象
    
    返回:
        SBtabDocument 对象
    """
```

---

## 常量与默认值

### 9. 包级常量

```python
# 默认相态
default_phase = "aqueous"

# 默认生理条件
default_physiological_p_h         # 默认生理 pH
default_physiological_p_mg        # 默认生理 pMg
default_physiological_ionic_strength  # 默认生理离子强度
default_physiological_temperature     # 默认生理温度

# 默认浓度边界
default_conc_lb   # 默认浓度下界
default_conc_ub   # 默认浓度上界

# 其他默认值
default_e_potential  # 默认电势
default_rmse_inf     # 默认 RMSE 无穷大值
```

---

## 完整使用示例

### 示例 1: 基本反应计算

```python
from equilibrator_api import ComponentContribution

# 初始化
cc = ComponentContribution()

# 解析反应
rxn = cc.parse_reaction_formula("ATP + H2O => ADP + Pi")

# 计算标准吉布斯自由能
dg0_prime = cc.standard_dg_prime(rxn)
print(f"ΔG'° = {dg0_prime.value:.2f} ± {dg0_prime.error:.2f} kJ/mol")

# 计算生理条件下的吉布斯自由能
dgm_prime = cc.physiological_dg_prime(rxn)
print(f"ΔG'm = {dgm_prime.value:.2f} ± {dgm_prime.error:.2f} kJ/mol")
```

### 示例 2: 设置自定义浓度

```python
from equilibrator_api import ComponentContribution, Q_

cc = ComponentContribution()
rxn = cc.parse_reaction_formula("glucose + ATP => glucose-6-phosphate + ADP")

# 获取化合物
atp = cc.get_compound("KEGG:C00002")
adp = cc.get_compound("KEGG:C00008")

# 设置自定义浓度
rxn.set_abundance(atp, 5 * Q_("mM"))
rxn.set_abundance(adp, 0.5 * Q_("mM"))

# 计算调整后的 ΔG'
dg_prime = cc.dg_prime(rxn)
print(f"ΔG' = {dg_prime.value:.2f} ± {dg_prime.error:.2f} kJ/mol")
```

### 示例 3: 批量计算多个反应

```python
from equilibrator_api import ComponentContribution

cc = ComponentContribution()

# 定义多个反应
formulas = [
    "glucose => 2 pyruvate",
    "pyruvate + CoA + NAD => acetyl-CoA + CO2 + NADH",
    "acetyl-CoA + 2 H2O + 3 NAD + FAD => 2 CO2 + 3 NADH + FADH2 + CoA"
]

# 解析反应
reactions = [cc.parse_reaction_formula(f) for f in formulas]

# 批量计算
dg_primes, uncertainties = cc.standard_dg_prime_multi(
    reactions,
    uncertainty_representation='cov'
)

# 显示结果
for i, (formula, dg) in enumerate(zip(formulas, dg_primes)):
    print(f"{formula}")
    print(f"  ΔG'° = {dg:.2f} kJ/mol")
    print(f"  Uncertainty: ±{uncertainties[i,i]**0.5:.2f} kJ/mol\n")
```

### 示例 4: 使用 Bounds 设置浓度范围

```python
from equilibrator_api import ComponentContribution, Q_
from equilibrator_api.model import Bounds

cc = ComponentContribution()

# 创建自定义边界
bounds = Bounds(
    default_lb=1e-6 * Q_("M"),
    default_ub=1e-2 * Q_("M")
)

# 为特定化合物设置边界
atp = cc.get_compound("KEGG:C00002")
bounds.set_bounds(atp, lb=1e-3 * Q_("M"), ub=5e-3 * Q_("M"))

# 获取边界
lb, ub = bounds.get_bounds([atp])
print(f"ATP concentration range: {lb[0]} to {ub[0]}")
```

### 示例 5: 电化学计算

```python
from equilibrator_api import ComponentContribution

cc = ComponentContribution()

# 定义半反应（还原反应）
half_rxn = cc.parse_reaction_formula("NAD + 2 e- + H+ => NADH")

# 计算标准还原电势
e0_prime = cc.standard_e_prime(half_rxn)
print(f"E'° = {e0_prime.value:.3f} ± {e0_prime.error:.3f} V")

# 计算生理条件下的还原电势
em_prime = cc.physiological_e_prime(half_rxn)
print(f"E'm = {em_prime.value:.3f} ± {em_prime.error:.3f} V")
```

### 示例 6: pH 敏感性分析

```python
from equilibrator_api import ComponentContribution
import numpy as np
import matplotlib.pyplot as plt

cc = ComponentContribution()
rxn = cc.parse_reaction_formula("ATP + H2O => ADP + Pi")

# 计算 pH 敏感性
sensitivity = cc.dg_prime_sensitivity_to_p_h(rxn)
print(f"∂(ΔG')/∂pH = {sensitivity:.2f} kJ/mol per pH unit")

# 绘制 ΔG' vs pH 曲线
ph_range = np.linspace(5, 9, 50)
dg_values = []

for ph in ph_range:
    cc.p_h = ph
    dg = cc.standard_dg_prime(rxn)
    dg_values.append(dg.value.magnitude)

plt.plot(ph_range, dg_values)
plt.xlabel('pH')
plt.ylabel("ΔG'° (kJ/mol)")
plt.title('ATP Hydrolysis: ΔG\'° vs pH')
plt.grid(True)
plt.show()
```

### 示例 7: 从 SBtab 文件加载模型

```python
from equilibrator_api import ComponentContribution
from equilibrator_api.model import StoichiometricModel

cc = ComponentContribution()

# 从 SBtab 文件加载
model = StoichiometricModel.from_sbtab(
    "pathway.tsv",
    comp_contrib=cc
)

# 查看化合物
print("Compounds:")
for cpd_id in model.compound_ids:
    print(f"  {cpd_id}")

# 查看反应
print("\nReactions:")
for rxn_id, formula in zip(model.reaction_ids, model.reaction_formulas):
    print(f"  {rxn_id}: {formula}")

# 获取反应数据框
df = model.reaction_df
print("\n", df)
```

### 示例 8: COBRA 模型集成

```python
from equilibrator_api import ComponentContribution
from equilibrator_api.compatibility import map_cobra_reactions
import cobra

# 加载 COBRA 模型
cobra_model = cobra.io.load_json_model("e_coli_core.json")

# 初始化 eQuilibrator
cc = ComponentContribution()

# 映射反应
equilibrator_reactions = map_cobra_reactions(
    cache=cc.ccache,
    reactions=cobra_model.reactions
)

# 计算热力学
for rxn_id, eq_rxn in equilibrator_reactions.items():
    try:
        dg_prime = cc.standard_dg_prime(eq_rxn)
        print(f"{rxn_id}: ΔG'° = {dg_prime.value:.2f} kJ/mol")
    except Exception as e:
        print(f"{rxn_id}: Error - {e}")
```

---

## 错误处理

### 常见错误及解决方案

#### 1. 化合物未找到

```python
from equilibrator_api import ComponentContribution

cc = ComponentContribution()

try:
    cpd = cc.search_compound("unknown_compound")
    if cpd is None:
        print("Compound not found")
except Exception as e:
    print(f"Error: {e}")
```

#### 2. 反应不平衡

```python
from equilibrator_api import ComponentContribution

cc = ComponentContribution()

try:
    rxn = cc.parse_reaction_formula("glucose => pyruvate")  # 不平衡
    # 使用氧化平衡
    balanced_rxn = cc.balance_by_oxidation(rxn)
    print(f"Balanced: {balanced_rxn}")
except Exception as e:
    print(f"Error: {e}")
```

#### 3. 边界检查

```python
from equilibrator_api.model import Bounds
from equilibrator_api import Q_

bounds = Bounds()
bounds.set_bounds("KEGG:C00002", lb=5 * Q_("mM"), ub=1 * Q_("mM"))  # 错误：lb > ub

try:
    bounds.check_bounds()
except AssertionError:
    print("Invalid bounds: lower bound must be <= upper bound")
```

---

## 高级主题

### 1. 自定义不确定性表示

```python
from equilibrator_api import ComponentContribution

cc = ComponentContribution()
reactions = [
    cc.parse_reaction_formula("glucose + ATP => glucose-6-phosphate + ADP"),
    cc.parse_reaction_formula("fructose-6-phosphate => fructose-1,6-bisphosphate")
]

# 使用不同的不确定性表示
dg_cov, cov_matrix = cc.standard_dg_prime_multi(reactions, uncertainty_representation='cov')
dg_sqrt, sqrt_matrix = cc.standard_dg_prime_multi(reactions, uncertainty_representation='sqrt')
dg_full, full_matrix = cc.standard_dg_prime_multi(reactions, uncertainty_representation='fullrank')

print("Covariance matrix shape:", cov_matrix.shape)
print("Square root matrix shape:", sqrt_matrix.shape)
print("Full-rank matrix shape:", full_matrix.shape)
```

### 2. 多室反应

```python
from equilibrator_api import ComponentContribution, Q_

cc = ComponentContribution()

# 定义内室和外室半反应
rxn_inner = cc.parse_reaction_formula("4 H+[in] => 4 H+[out]")
rxn_outer = cc.parse_reaction_formula("ATP + H2O => ADP + Pi")

# 计算跨膜转运的 ΔG'
dg_transport = cc.multicompartmental_standard_dg_prime(
    reaction_inner=rxn_inner,
    reaction_outer=rxn_outer,
    e_potential_difference=150 * Q_("mV"),
    p_h_outer=7.0,
    ionic_strength_outer=0.25 * Q_("M")
)

print(f"Transport ΔG' = {dg_transport.value:.2f} kJ/mol")
```

### 3. 使用旧版本（向后兼容）

```python
from equilibrator_api import ComponentContribution

# 使用旧版本以保持与旧代码的兼容性
cc_legacy = ComponentContribution.legacy()

rxn = cc_legacy.parse_reaction_formula("ATP + H2O => ADP + Pi")
dg_legacy = cc_legacy.standard_dg_prime(rxn)

print(f"Legacy ΔG'° = {dg_legacy.value:.2f} kJ/mol")
```

---

## 性能优化建议

### 1. 批量计算优先

```python
# ❌ 不推荐：逐个计算
for rxn in reactions:
    dg = cc.standard_dg_prime(rxn)

# ✅ 推荐：批量计算
dg_values, uncertainties = cc.standard_dg_prime_multi(reactions)
```

### 2. 重用 ComponentContribution 对象

```python
# ✅ 推荐：重用对象
cc = ComponentContribution()
for formula in formulas:
    rxn = cc.parse_reaction_formula(formula)
    dg = cc.standard_dg_prime(rxn)
```

### 3. 缓存化合物查询

```python
# ✅ 推荐：缓存常用化合物
compound_cache = {}
for cpd_id in ["KEGG:C00002", "KEGG:C00008", "KEGG:C00020"]:
    compound_cache[cpd_id] = cc.get_compound(cpd_id)
```

---