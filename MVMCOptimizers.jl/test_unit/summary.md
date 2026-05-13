# test_unit 作業ログ

- テスト一覧: `INDEX.md` を参照（`test_unit/INDEX.md`）

- 日時: 2026-02-06 17:56:40 JST
- AIモデル: GPT-5 (Codex)

## 実施内容
1. `test_unit` 方針書を作成。
- `test_unit/plan.md` を追加し、対象関数・対象ファイル・優先順位を整理。

2. 1ファイル目の unit test を実装。
- `test_unit/test_unit_stochastic_opt.jl` を新規作成。
- 対象関数:
- `get_opt_flag_for_parameter`
- `update_parameter_value`
- `build_s_matrix_and_g_vector!`
- 主な検証:
- `Proj / RBM / Slater` のパラメータ更新マップ
- RBM 9セクション更新の整合
- `S` / `g` の手計算一致

3. テスト実行導線を追加。
- `test/runtests.jl` に `@testset "Unit Tests"` を追加。
- `include("../test_unit/test_unit_stochastic_opt.jl")` を追加。

## 実行結果
- 単体実行:
- `julia --project=@. -e 'include("test_unit/test_unit_stochastic_opt.jl")'`
- 結果: 全PASS（`6 + 13 + 6` tests）

- 全体実行:
- `julia --project=@. -e 'include("test/runtests.jl")'`
- 結果: `Unit Tests` は `25/25 PASS`、全体も終了コード `0` で完走

## 変更ファイル
- `test_unit/plan.md`
- `test_unit/test_unit_stochastic_opt.jl`
- `test/runtests.jl`

## 再開メモ
- 単体確認:
- `julia --project=@. -e 'include("test_unit/test_unit_stochastic_opt.jl")'`
- 全体確認:
- `julia --project=@. -e 'include("test/runtests.jl")'`
- 次の着手ファイル:
- `test_unit/test_unit_vmc_sampling_rbm.jl`

---

- 日時: 2026-02-14
- AIモデル: GPT-5.2 (Codex CLI)

## 実施内容（追記）
1. unit test を追加・拡充。
- `test_unit/test_unit_vmc_sampling_rbm.jl`（RBM カウンタの incremental update / full recompute の一致、`log_rbm_ratio` 整合）
- `test_unit/test_unit_vmc_sampling_proj.jl`（proj カウンタの incremental update / full recompute の一致、`fsz` variant も確認）
- `test_unit/test_unit_vmc_main_cal_sr.jl`（`set_projection_diff!` / `set_rbm_diff!` / `calculate_oo!` / `calculate_oo_real!`）
- `test_unit/test_unit_parameter_sync.jl`（`sync_modified_parameter!` の shift / rescale / normalize）

2. helper を導入し、fixture と move を集約。
- `test_unit/helpers/mock_state.jl`（`make_ele_num` / `apply_hop`）
- `test_unit/helpers/mock_data.jl`（`ExpertModeData` の最小 fixture）

3. unit test 実行導線を更新。
- `test/runtests.jl` の `Unit Tests` 内で helper を include してから各 `test_unit/*.jl` を include。

## 実行結果
- `julia --project=@. -e 'import Pkg; Pkg.test()'`
- 結果: `Unit Tests 118/118 PASS`、`Pkg.test()` 全体も終了コード `0` で完走（integration はデフォルト無効で `@test_skip`）。
