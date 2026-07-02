# multi_source_compare.jl
# 比较三组概率分布下三种编码方法的性能

include("source_coding_functions.jl")
using Printf

p_set = [
    [0.25, 0.25, 0.25, 0.25],
    [0.40, 0.30, 0.20, 0.10],
    [0.30, 0.25, 0.20, 0.15, 0.10]
]

count_set = [
    [5, 5, 5, 5],
    [8, 6, 4, 2],
    [6, 5, 4, 3, 2]
]

println("不同概率分布下的编码性能比较")
println("编号\t信源熵\tHuffman平均码长\tHuffman效率\tShannon平均码长\tShannon效率\t算术平均码长\t算术效率")

for case_index in eachindex(p_set)
    p = p_set[case_index]
    counts = count_set[case_index]
    symbols = collect(1:length(p))

    source = Int[]
    for i in eachindex(counts)
        append!(source, fill(symbols[i], counts[i]))
    end

    H = -sum(p .* log2.(p))

    _, huffman_lengths = huffman_build(p)
    huffman_average_length = sum(p .* huffman_lengths)
    huffman_efficiency = H / huffman_average_length

    _, shannon_lengths = shannon_build(p)
    shannon_average_length = sum(p .* shannon_lengths)
    shannon_efficiency = H / shannon_average_length

    arithmetic_bits, _, _, _ = arithmetic_encode(source, symbols, p)
    arithmetic_average_length = length(arithmetic_bits) / length(source)
    arithmetic_efficiency = H / arithmetic_average_length

    @printf("%d\t%.4f\t%.4f\t\t\t%.4f%%\t\t%.4f\t\t\t%.4f%%\t\t%.4f\t\t%.4f%%\n",
            case_index,
            H,
            huffman_average_length,
            huffman_efficiency * 100,
            shannon_average_length,
            shannon_efficiency * 100,
            arithmetic_average_length,
            arithmetic_efficiency * 100)
end
