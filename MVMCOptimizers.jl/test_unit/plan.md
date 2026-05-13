# Unit Test Plan for J-mVMC

## 目的
- C/J 比較の統合テストで見つかった差分を、関数レベルで早期検知できるようにする。
- `NProj + NRBM + NSlater` のようなインデックス規約や、fsz/cmp 分岐の回帰を unit test で固定する。

## テスト階層
- unit: 小関数・純関数・小規模状態（外部ファイルなし、1秒以内）
- contract: C 実装との対応仕様（小さな fixture を使用）
- integration: 既存の `test/test_vmc_models.jl`（モデル全体比較）

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
- `test/runtests.jl`:
- unit test は常時実行
- 統合テスト (`test_vmc_models.jl`) は `MVMC_INTEGRATION_TESTS=1` のときのみ実行

## 実装順序（優先度）
1. `test_unit_stochastic_opt.jl`
2. `test_unit_vmc_sampling_rbm.jl`
3. `test_trans_parser_spin_indices.jl`
4. `test_read_input_parameters_rbm_layout.jl`
5. `test_unit_vmc_main_cal_sr.jl`

## 再開コマンド
- 単体確認:
- `julia --project=@. -e 'include("test_unit/test_unit_stochastic_opt.jl")'`
- 全体確認:
- `julia --project=@. -e 'include("test/runtests.jl")'`

## 次の着手ファイル
- （要相談）次の unit/contract テスト対象を選定
