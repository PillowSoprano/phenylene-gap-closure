# 作图数据（Figs. 1–3）

所有图只从这些 CSV 画出：`plot_figs.py` 只读这里，不做任何计算。数据由 `export_figdata.py` 生成，计算与 `make_figs.py` 完全相同；Fig. 3 直接取自 `dft_validation/*_points.json`。想用 TikZ/pgfplots 重画，直接 `\addplot table[col sep=comma, x=..., y=...]{fig_data/xxx.csv}` 即可。

坐标单位：分子坐标以键长为 1（示意坐标，不是优化几何）。能量单位 eV，能隙在 Fig. 1 中以 |β| 为单位。

| 文件 | 列 | 内容 |
|---|---|---|
| fig1_atoms.csv | atom, x, y | 联苯撑原子坐标 |
| fig1_bonds.csv | atom1, atom2, cooperative_sign, single_bond_route | 每条键；cooperative_sign = +1 加强 / −1 削弱（阈值处）；single_bond_route = 1 表示单键路线要放大的那根键 |
| fig1_curves.csv | route, largest_relative_change_percent, huckel_gap_abs_beta | 两条路线（single / cooperative）上的 Hückel 能隙 |
| fig1_summary.csv | quantity, value | ρ(联苯撑)，单键路线所需的相对增量 |
| fig2_<skeleton>_atoms.csv | atom, x, y, closing_orbital_weight | 闭合轨道在每个原子上的权重（全分子求和为 1） |
| fig2_<skeleton>_bonds.csv | atom1, atom2, coupling_change | 阈值处每条键耦合的相对变化（±ρ） |
| fig2_summary.csv | skeleton, rho, sigma_min_at_pattern | 精确阈值；闭合模式处的最小奇异值（应为机器零） |
| fig3_ladders.csv | molecule, pattern, model_energy_eV, b3lyp_deformation_energy_eV, b3lyp_gap_eV, gap_decrease_eV, y0, S2_broken_symmetry | 预测模式（P）与键长均匀化对照（BLA）的阶梯；参考点在两条阶梯中各出现一次；y0 和 ⟨S²⟩ 只在 P 阶梯上有值 |
| fig3_random.csv | molecule, label, b3lyp_deformation_energy_eV, b3lyp_gap_eV, gap_decrease_eV | 16 个随机模式（8 个方向 × 正负号），模型能量均为 1.5 eV |
