# Unit Test Plan for J-mVMC

## 目的
- C/J 比較の統合テストで見つかった差分を、関数レベルで早期検知できるようにする。
- `NProj + NRBM + NSlater` のようなインデックス規約や、fsz/cmp 分岐の回帰を unit test で固定する。

## テスト階層
- unit: 小関数・純関数・小規模状態（外部ファイルなし、1秒以内）
- contract: 入力仕様・runtime 互換性の固定（小さな fixture / 入力検証）
- integration: C reference との突き合わせ。public repo では workspace root の
  `test/integration/runtests.jl` が担当（`julia --project=@. test/integration/runtests.jl`）。
  subpackage の `Pkg.test()` には含まれない別レイヤ。

## 追加対象ファイル（MVMCOptimizers.jl）

### 1. `test_unit/test_unit_stochastic_opt.jl`
- 対象: `src/stochastic_opt.jl`
- 関数:
- `get_opt_flag_for_parameter`
- `update_parameter_value`
- `build_s_matrix_and_g_vector!`
- 主な検証:
- `update_parameter_value` の `Proj / RBM / Slater` マップが正しいこと
- RBM 9セクション（phys(3) -> hidden(3) -> phys-hidden(3)）更新が正しいこと
- `build_s_matrix_and_g_vector!` の `S` と `g` が期待式どおりであること

### 2. `test_unit/test_unit_vmc_sampling_rbm.jl`
- 対象: `src/vmc_sampling.jl`
- 関数:
- `has_rbm_terms`
- `make_rbm_cnt`
- `update_rbm_cnt_hopping!`
- `log_rbm_ratio`
- `log_rbm_val`
- 主な検証:
- `incremental update` と `full recompute` の一致
- 1 hop / 2 hop (exchange) でカウンタ更新が一致
- `log_rbm_ratio(new, old)` が `log_rbm_val(new_ele) - log_rbm_val(old_ele)` と整合

### 3. `test_unit/test_unit_vmc_sampling_proj.jl`
- 対象: `src/vmc_sampling.jl`
- 関数:
- `make_proj_cnt!`
- `update_proj_cnt!`
- `update_proj_cnt_fsz!`
- 主な検証:
- update 前後で `make_proj_cnt!` 再計算結果と一致
- fsz と非fszで同じ更新規約を満たすこと

### 4. `test_unit/test_unit_vmc_main_cal_sr.jl`
- 対象: `src/vmc_main_cal.jl`
- 関数:
- `set_projection_diff!`
- `set_rbm_diff!`
- `calculate_oo!`
- `calculate_oo_real!`
- 主な検証:
- SR ベクトル配置（オフセット）が仕様どおり
- `OO/HO` 蓄積で実装式と期待値が一致

### 5. `test_unit/test_unit_parameter_sync.jl`
- 対象: `src/parameter_sync.jl`
- 関数:
- `sync_modified_parameter!`
- 主な検証:
- Slater の最大振幅が `D_AMP_MAX` に収まること
- Gutzwiller/Jastrow shift 条件が仕様どおりであること

## 追加対象ファイル（MVMCExpertModeParsers.jl）

### 6. `test/test_trans_parser_spin_indices.jl`
- 対象: `src/parsers/trans_parser.jl`
- 関数:
- `parse_trans_content`
- `parse_transfer_term`
- 主な検証:
- `spin1/spin2` が正しく保持されること
- `spin1 != spin2` ケースの扱いが仕様どおりであること

### 7. `test/test_read_input_parameters_rbm_layout.jl`
- 対象: `src/utils/read_input_parameters.jl`
- 関数:
- `count_rbm_parameters`
- `set_rbm_opt_flags!`
- `read_input_parameters!`
- 主な検証:
- RBM 9セクション順の offset が崩れないこと
- OptFlag サイズと割り当てが正しいこと

### 8. `test/test_parameter_init_complexflag_rbm.jl`
- 対象: `src/utils/parameter_init.jl`
- 関数:
- `init_parameter!`
- `initialize_parameters!`
- 主な検証:
- `AllComplexFlag` 相当判定が期待どおり
- RNG seed 固定時に初期化順と結果が再現すること

## テストデータ設計
- 新規 helper:
- `test_unit/helpers/mock_data.jl`
- `test_unit/helpers/mock_state.jl`
- 用意する生成関数:
- `make_minimal_data_real()`
- `make_minimal_data_fsz()`
- `make_minimal_data_rbm()`
- `make_minimal_state(data; n_sample=...)`

## 実行導線
- subpackage `test/runtests.jl`:
- unit test（`test_unit/*.jl`）を常時実行。include 順は `INDEX.md` と同期する。
- C reference との integration は subpackage には含めない。workspace root の
  `test/integration/runtests.jl` で別途実行する（`julia --project=@. test/integration/runtests.jl`）。
- 対応 Julia は 1.11+（`Project.toml` compat `julia = "1.11"`、CI は 1.11 / 1.12 を検証）。
  `Manifest.toml` は gitignore 対象。ローカルで 1.11 / 1.12 を切り替える場合は、その version で
  resolve（`Pkg.instantiate` / `Pkg.resolve`）し直すこと。

## 状態（v0.1 計画ぶんの達成状況）
- 上記「追加対象ファイル」の unit/contract テストは全て実装・常時実行済み:
  `test_unit_stochastic_opt.jl` / `test_unit_vmc_sampling_rbm.jl` /
  `test_unit_vmc_sampling_proj.jl` / `test_unit_vmc_main_cal_sr.jl` /
  `test_unit_parameter_sync.jl`、および parsers 側
  `test_trans_parser_spin_indices.jl` / `test_read_input_parameters_rbm_layout.jl` /
  `test_parameter_init_complexflag_rbm.jl`。
- helper は `helpers/mock_state.jl`（`make_ele_num` / `apply_hop`）と
  `helpers/mock_data.jl`（最小 `ExpertModeData` constructors）に集約済み。

## v0.3 優先度
詳細設計は `docs/superpowers/specs/2026-06-03-julia-mvmc-v0.3-unit-test-foundation-design.md`。

1. [済] NSplitSize contract: runtime に `validate_supported_modpara` を追加し
   `NSplitSize > 1` を `error()`/`ErrorException` で reject。
   test: `test_unit/test_unit_unsupported_inputs.jl`。
2. PhysCal production 対応 + PhysCal C-reference fixtures。
3. DH2/DH4・InDH2/InDH4 と warn-only 入力互換。
4. OpenMP 相当の Julia threading（thread-local accumulator / timer / RNG state）。

## 再開コマンド
- 単体確認:
- `julia --project=@. -e 'include("test_unit/<file>.jl")'`
- 全体確認:
- `julia --project=@. -e 'using Pkg; Pkg.test()'`
