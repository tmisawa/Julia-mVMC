# Unit Test Index (MVMCOptimizers.jl)

このフォルダ配下の unit test は `test/runtests.jl` から常時実行されます。C reference との integration は subpackage には含めず、workspace root の `test/integration/runtests.jl` で別途実行します。対応 Julia は 1.11+（`Project.toml` compat `julia = "1.11"`、CI は 1.11 / 1.12 を検証）。`Manifest.toml` は gitignore 対象なので、ローカル実行時は使用する Julia version で resolve し直すこと。

## 実行
- 全体: `cd MVMCOptimizers.jl && julia --project=@. -e 'import Pkg; Pkg.test()'`（リポジトリルートから）
- unit のみ（テストファイル単体）: `julia --project=@. -e 'include("test_unit/<file>.jl")'`

## 読み込み経路
- `test/runtests.jl` の `@testset "Unit Tests"` で include（`test/runtests.jl:62` 付近）

## Helper（fixture / state）
- `test_unit/helpers/mock_state.jl`
  - `make_ele_num(nsite; up_sites, down_sites)`（`ele_num` 生成）
  - `apply_hop(ele_num, ri, rj, s, nsite)`（単一 hop 適用）
- `test_unit/helpers/mock_data.jl`
  - `make_mock_data_for_rbm_tests(; nsite)`（RBM 用の最小 `ExpertModeData`）
  - `make_mock_data_for_proj_tests(; nsite)`（Proj 用の最小 `ExpertModeData`）
  - `make_minimal_data_for_rbm_diff_tests(; nsite, nneuron_charge)`（SR diff 用の最小 `ExpertModeData`）

## Unit Tests 一覧

### `src/stochastic_opt.jl`
- `test_unit/test_unit_stochastic_opt.jl`
  - `unit/stochastic_opt: update_parameter_value DH write-back and shifted offsets` → `update_parameter_value`（DH2 projection block と DH 後続 RBM/Slater offset）
  - `unit/stochastic_opt: update_parameter_value DH4 write-back` → `update_parameter_value`（DH4 projection block と DH4 後続 RBM/Slater offset）
  - `unit/stochastic_opt: SR enumeration writes DH without touching RBM or Slater` → `stochastic_opt!` / `update_parameter_value`
  - `unit/stochastic_opt: get_opt_flag_for_parameter` → `get_opt_flag_for_parameter`
  - `unit/stochastic_opt: update_parameter_value Proj/RBM/Slater mapping` → `update_parameter_value`
  - `unit/stochastic_opt: build_s_matrix_and_g_vector!` → `build_s_matrix_and_g_vector!`

### `src/vmc_sampling.jl`
- `test_unit/test_unit_vmc_sampling_rbm.jl`
  - `unit/vmc_sampling: has_rbm_terms` → `has_rbm_terms`
  - `unit/vmc_sampling: make_rbm_cnt / update_rbm_cnt_hopping! / log_rbm_* consistency`
    - `make_rbm_cnt`
    - `update_rbm_cnt_hopping!`
    - `log_rbm_ratio`
    - `log_rbm_val`

- `test_unit/test_unit_vmc_sampling_proj.jl`
  - `unit/vmc_sampling: make_proj_cnt! / update_proj_cnt! consistency`
    - `make_proj_cnt!`
    - `update_proj_cnt!`
    - `update_proj_cnt_fsz!`
  - DH projection count / update / log coverage
    - `unit/vmc_sampling: DH2/DH4 projection counts` → `make_proj_cnt!` / `set_projection_diff!`
    - `unit/vmc_sampling: multiple DH2 indices keep independent strides` → `make_proj_cnt!`（`n_dh2 >= 2` stride）
    - `unit/vmc_sampling: DH update paths match fresh recompute` → `update_proj_cnt!` / `update_proj_cnt_fsz!`（DH2, alias buffer, double-hop alias）
    - `unit/vmc_sampling: DH4 update paths match fresh recompute` → `update_proj_cnt!` / `update_proj_cnt_fsz!`（DH4, alias buffer）
    - `unit/vmc_sampling: DH log projection uses real parameter parts` → `log_proj_val` / `log_proj_ratio`

- `test_unit/test_unit_vmc_sampling_misc.jl`
  - `unit/vmc_sampling: log_proj_val / log_proj_ratio`
    - `log_proj_val`
    - `log_proj_ratio`
  - `unit/vmc_sampling: update_ele_config! / revert_ele_config! round-trip`
    - `update_ele_config!`
    - `revert_ele_config!`

### `src/slater_update.jl`
- `test_unit/test_unit_slater_update.jl`
  - `unit/slater_update: build_orbital_idx_sgn_matrices` → `build_orbital_idx_sgn_matrices`
  - `unit/slater_update: build_orbital_idx_sgn_matrices uses data.orbital_sgn when provided` → 同上（`data.orbital_sgn` 優先）
  - `unit/slater_update: build_qp_trans_matrices requires qp_trans mappings` → `build_qp_trans_matrices`（必須マップ欠落時の挙動）
  - `unit/slater_update: build_qp_trans_matrices opt fallback is identity` → `build_qp_trans_matrices`（`qp_opt_trans*` の fallback）

### `src/vmc_main_cal.jl`
- `test_unit/test_unit_vmc_main_cal_sr.jl`
  - `unit/vmc_main_cal: set_projection_diff!` → `set_projection_diff!`
  - `unit/vmc_main_cal: set_rbm_diff!` → `set_rbm_diff!`
  - `unit/vmc_main_cal: calculate_oo!` → `calculate_oo!`
  - `unit/vmc_main_cal: calculate_oo_real!` → `calculate_oo_real!`
  - `unit/vmc_main_cal: threaded OO/store inner loops match sequential` → R3 inner-loop threading（complex/real OO, store, finalize）

### `src/parameter_sync.jl`
- `test_unit/test_unit_parameter_sync.jl`
  - `unit/parameter_sync: sync_modified_parameter!` → `sync_modified_parameter!`
    - shift（Gutzwiller/Jastrow）
    - shift（DH2/DH4, Gutzwiller write-back, OptFlag guard）
    - rescale（`D_AMP_MAX`）
    - normalize（`opt_trans`、`para_qp_trans` 非正規化ガード）
    - shift（DH2+DH4 同時 config）

### `src/types.jl`
- `test_unit/test_unit_types.jl`
  - `unit/types: EnergyData defaults` → `EnergyData()`
  - `unit/types: SROptData allocation sizes` → `SROptData(sr_opt_size, n_vmc_sample, all_complex)`
  - `unit/types: ElectronConfiguration sizes (fsz vs non-fsz)` → `ElectronConfiguration(..., use_fsz)`
  - `unit/types: SlaterMatrixData size normalization` → `SlaterMatrixData(n_qp_full, n_site, n_elec, all_complex)`

### `src/threading.jl`（inner-loop threading helpers / local accumulators）
- `test_unit/test_unit_threading.jl`
  - `unit/threading: inner copy helpers` → R3 real/complex copy helper の threaded/serial 一致
  - `unit/threading: energy accumulator reduction` → deterministic local energy merge
  - `unit/threading: SR accumulator reduction` → complex/real SR array merge
  - `unit/threading: PhysCal accumulator reduction` → `phys_*` merge without touching local scratch、factored Green の local accumulator 加算
  - `unit/threading: counter and timer reductions` → counter sum and `CTimer` elapsed merge
  - `unit/threading: local accumulator construction and merge` → serial VMCMainCal local accumulator allocation and parent-state merge

### `src/unsupported_inputs.jl`（runtime contract）
- `test_unit/test_unit_unsupported_inputs.jl`
  - `unit/unsupported_inputs: NSplitSize contract` → `validate_supported_modpara`
    - `NSplitSize = 1` は許容（serial / v0.4 MPI sample parallel、`vmc_para_opt!` / `vmc_phys_cal!`）
    - `NSplitSize > 1` は `error()` で reject（C の grouped MPI/QP split 未サポート）
    - エラーメッセージに `NSplitSize > 1` と `grouped MPI/QP splitting by NSplitSize is not implemented` を含む
    - 検証は型ではなくメッセージ部分文字列で行う（design review A2）
  - `unit/unsupported_inputs: NLanczosMode contract` → `validate_supported_modpara`
    - `NLanczosMode = 0` は許容、`> 0`（1/2）は `error()` で reject（full Lanczos 未サポート / C の indirect one-body list 非再現）
    - `vmc_para_opt!` / `vmc_phys_cal!` の両 entry point で enforce、メッセージに `NLanczosMode > 0` を含む

### `src/green_func_calc.jl` + `vmc_phys_cal.jl`（factored two-body Green）
- `test_unit/test_unit_physcal_factored_green.jl`
  - `PhysicalQuantities index fields` → `cis_ajs_idx` / `cis_ajs_ckt_alt_idx`
  - `canonical one-body list` → `build_canonical_cis_ajs_idx`（greenone 先頭 → greentwoex 構成 append・dedup・site 範囲・spin 検証）
  - `factored index resolution is 1-based` → `resolve_cis_ajs_ckt_alt_idx`（C index 0 → Julia 1）
  - `factored accumulation` → `accumulate_factored_green!`（`w·local[idx0]·conj(local[idx1])`）
  - `calculate_green_func! legacy no-acc path writes to state phys` → accumulator なし旧 API の R1 follow-up 回帰ガード
  - `output: canonical cisajs + factored ex` → `output_green_func!`（canonical 出力 / output_dir / `_001` 番号）
  - `PhysCal output file index uses NDataIdxStart` → `physcal_output_file_index`
  - `output_data_phys!: out/var truncate-on-first-sample, Green files use NDataIdxStart` → fmt-1 回帰ガード（out/var は 0-based ismp で write-mode 決定 → 初回 truncate・再実行で非汚染、Green は `ismp+NDataIdxStart` 番号）
  - `initialize_phys_quantities! wires the factored canonical list and pairs` → factored 分岐の統合（canonical/pairs を struct へ配線・buffer サイズ、共有構成2項の dedup append + 複数ペア順序）
  - `no TwoBodyGEx preserves greenone order and duplicates` → 既存 direct 経路の互換回帰
  - `FSZ + factored is rejected` → `validate_factored_green_supported`（public path `vmc_phys_cal!` 経由も検証）

### `src/initial_params.jl`（fixed-parameter loader, Plan 3a）
- `test_unit/test_unit_read_opt_para.jl`
  - `read_opt_para_file!: golden load + consumed count` → 6+triples レイアウトを Gutzwiller/Jastrow/Slater（no-DH）に読み込み、消費数（n_proj+n_slater）を返す
  - `non-perturbing (idempotent) load` → 二重ロードで同一（RNG 非依存）
  - `strict failures error with a clear message` → 欠損 / 空 / short / trailing（1 float・inactive OptTrans whole triple）/ non-numeric token を `error()`（メッセージ部分文字列で検証）
  - `DH projection block loads before Slater` → `NProj` 内の DH2 block を読み込み、Slater が DH の後ろから正しく scatter されることを検証
  - `OptTrans tail loads when active` → `opttrans.def` active 相当の `data.opt_trans` がある場合、Slater 後ろの `NOptTrans` triple を消費
  - `RBM-bearing models still fail loud` → 将来 RBM block を Slater と誤読しないため、`NRBM > 0` は DH-2 でも loader reject
  - `read_initial_def! still loads the same layout (delegation regression)` → 共有 `_load_para_triples!` への委譲後も成功/欠損 warn+false が不変

### `src/run_phys_cal_from_namelist.jl`（PhysCal runner, Plan 3a）
- `test_unit/test_unit_run_phys_cal_runner.jl`（Pkg.test）
  - `argument validation` → mode（:real/:cmp/:fsz）検証・`opt_para` は必須 kwarg。runner は `init_parameter!`/`init_qp_weight!` を呼ばず `vmc_phys_cal!` に委譲（RNG 二重消費回避）
- `../../test/integration/test_run_phys_cal_contract.jl`（CI: integration runtests 経由）
  - `errors at the loader, before sampling` → 既存 fixture(heisenberg_chain_real) で parse 成功後、欠損/short な `opt_para` を sampling 前に `read_opt_para_file!` error（出力 dir 未作成も確認）。数値 e2e は Plan 3b gate

### `test/integration/tools/green_compare.jl`（Green ファイル比較 helper, Plan 3a / 3b gate 用）
- `test/integration/tools/test_green_compare.jl`（CI: integration runtests 経由 + 単体実行可）
  - raw-byte exact → 数値 fallback（per-quantity: one-body rtol 1e-10 / factored・DC rtol 1e-9, atol 1e-12）。index 列は厳密一致、行/列/値数不一致・tol 超過は hard fail。raw-byte は trailing 空白に敏感（strip しない）。factored は単一行フォーマットを fallback 前に強制（layout regression guard）

### `test/integration/phys_cal_equivalent.jl`（PhysCal e2e gate, Plan 3b）
- CI: 専用 step `julia --project=@. test/integration/phys_cal_equivalent.jl`（`JULIA_MVMC_PHYS_CAL_MODELS` で subset 可）。subpackage の Pkg.test には含めない（explicit file invocation, ctest_equivalent.jl と同方式）
  - 5 系（heisenberg_chain_real/cmp, hubbard_chain_real, hubbard_chain_dh_real, kondo_chain_real）の `reference/<sys>/physcal_ref/`（NVMCCalMode=1 + 手書き `greentwoex.def`、DH 系は `dh2.def`/`dh4.def` も含む）に対し `run_phys_cal_from_namelist` を tempdir で実行 → `status==0`・`n_para_consumed` 厳密一致・3 Green 出力（one-body/direct/factored）を committed C reference と green_compare（`test_run_phys_cal_contract.jl` が defer した数値 e2e の実体）
  - C reference は主に mVMC develop @ 66f1742 を Apple Clang + gfortran で build（USE_GEMMT=OFF）。DH fixture は同 commit ベースの local benchmark branch @ 622166a で生成

---

## 参考: MVMCExpertModeParsers.jl 側の contract テスト
（MVMCOptimizers の unit test とは別パッケージですが、入力仕様の固定として関連します）
- `../MVMCExpertModeParsers.jl/test/test_trans_parser_spin_indices.jl`
  - `contract/trans.def: spin1/spin2 are preserved`
- `../MVMCExpertModeParsers.jl/test/test_read_input_parameters_rbm_layout.jl`
  - `contract/read_input_parameters: RBM layout (count/offsets)`
