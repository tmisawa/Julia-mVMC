# 概要

[English](../en/index.md)

Julia-mVMC は、多変数変分モンテカルロソルバ mVMC の Julia 実装です。
C-mVMC の expert-mode `.def` 入力を読み込み、C reference 実装との再現性を
重視しています。

この v0.5.0 manual は、現在の release で公開・検証済みの範囲を説明します。
C-mVMC の全機能を置き換えるものではありません。未対応の組合せは明示し、
静かに非互換な結果を出すのではなく、早い段階で error にする方針です。

## v0.5.0 の状態

| Component | Status | Notes |
|-----------|--------|-------|
| `VMCParaOpt` | 検証済み | first-10-step C-reference gate、C ctest-equivalent gate、serial / MPI `NSRCG = 1` first-step tolerance gate で確認済み。 |
| `VMCPhysCal` | 対応範囲は検証済み | one-body、direct two-body、factored/product two-body Green functions を対応 fixture で C reference と比較済み。 |
| Full Lanczos PhysCal | 対応範囲は検証済み | sz-conserved `NSplitSize = 1` path の `NLanczosMode = 1/2` output を C reference と比較済み。 |
| MPI parallelization | 部分検証 / smoke-tested | direct-SR split、standard-projection split、SR-CG `NSplitSize = 1`、rank0 output、PhysCal split を targeted smoke gate で確認。 |
| Shared-memory threading | experimental opt-in | conservative な inner-loop opt-in のみ。C parity のため sample-level Markov-chain threading は無効。 |
| BackFlow | 未対応 | 将来 release で対応予定。 |

## 波動関数の規約

Julia-mVMC は C-mVMC と同じ高レベルの変分波動関数を使います。

```math
|\Psi\rangle = \mathcal{P} |\Phi\rangle,
```

ここで ``|\Phi\rangle`` は Slater または generalized orbital の reference
state、``\mathcal{P}`` は Gutzwiller、Jastrow、doublon-holon、
quantum-number projection などの projection / correlation factor です。

sampling される変分エネルギーは次の形です。

```math
E = \frac{\langle \Psi | H | \Psi \rangle}
         {\langle \Psi | \Psi \rangle}.
```

v0.5.0 の実装目標は、新しい stochastic semantics を入れることではなく、
検証済み path で C reference と文書化された浮動小数点 tolerance 内で一致する
ことです。

## Manual の構成

- [インストール](installation.md): native 依存と clone-based setup。
- [チュートリアル](tutorial.md): bundled Heisenberg-chain optimization の実行。
- [入力ファイル](input_files.md): 対応 expert-mode `.def` file。
- [最適化](optimization.md): `VMCParaOpt` の使い方と MPI scope。
- [物理量計算](physics_calc.md): `VMCPhysCal` Green / Lanczos output。
- [出力ファイル](output_files.md): 主な `zvo_*` / `zqp_*` file。
- [互換性](compatibility.md): C-reference tests と未対応機能。
