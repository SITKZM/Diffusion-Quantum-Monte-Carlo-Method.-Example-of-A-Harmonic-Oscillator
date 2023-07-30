# Diffusion-Quantum-Monte-Carlo-Method.-Example-of-A-Harmonic-Oscillator
# 日本語
## このプログラムについて
このプログラムに関する詳細な解説は私のHPのノートを参照してください。

## 使い方
"DQMC.f90"に好きなパラメータを入れてコンパイルし、実行してください。
モンテカルロシミュレーションが始まり、終了すると3つのファイルが作成されます。
その3つのファイルの例を"Output_example"フォルダに入れておきました。"info.txt"には入力パラメータと出力データがまとめられており, "ER_evolution"にはシミュレーション前半における基準エネルギーの時間発展が、"phi_0.txt"にはシミュレーション結果として得られる基底状態波動関数のデータが入っています。後者2つの各行・列が何を表すかは、"DQMC.f90"内のサブルーチン"write_files()"を見れば分かると思います。

シミュレーションが終わったら、"plot.py"を実行することで結果の可視化が行えます。プロットされるのは基準エネルギー（＝基底状態エネルギー）の時間発展と基底状態波動関数です。どちらも灰色の線が引かれますが、それらは解析解です。

# English
## About this program
Please refer to the notes on my website for a detailed explanation of this program.

## Usage
Compile and run "DQMC.f90" with your desired parameters.
The Monte Carlo simulation will start, and when it finishes, three files will be created.
Examples of the three files are in the "Output_example" folder. "info.txt" contains a summary of input parameters and output data, "ER_evolution" contains the time evolution of the reference energy in the first half of the simulation, and "phi_0.txt" contains the resulting ground-state wave function. You can see what each row and column represents in the latter two files by looking at the subroutine "write_files()" in "DQMC.f90".

After the simulation, you can visualize the results by running "plot.py". Plotted are the time evolution of the reference energy (= ground state energy) and the ground state wave function. Both are plotted with gray lines, but they are analytical solutions.
