使用说明

1. MWORKS的数值计算环境采用Julia语法，脚本扩展名为.jl。
2. 将以下三个.jl文件放在同一目录：
   source_coding_functions.jl
   main_source_coding.jl
   multi_source_compare.jl
3. 首先运行main_source_coding.jl，获得主实验完整结果。
4. 再运行multi_source_compare.jl，获得不同概率分布的对比结果。
5. 程序只使用Julia基础语法和Printf标准库，不依赖额外工具箱。


