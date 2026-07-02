# main_source_coding.jl
# 主实验：同一离散无记忆信源采用三种编码方法

include("source_coding_functions.jl")
using Printf

symbols = [1, 2, 3, 4, 5]
p = [0.30, 0.25, 0.20, 0.15, 0.10]

source = [
    1, 2, 3, 4, 5,
    1, 2, 3, 1, 4,
    2, 1, 3, 5, 2,
    1, 4, 3, 2, 1
]

H = -sum(p .* log2.(p))
fixed_length = ceil(Int, log2(length(symbols)))
fixed_total_bits = fixed_length * length(source)

println("====================================================")
println("Huffman、Shannon和算术编码实验")
println("原始信源序列：", join(source, " "))
@printf("信源熵 H = %.6f bit/符号\n", H)
println("定长编码基准：", fixed_length, " bit/符号，共 ", fixed_total_bits, " bit")
println("====================================================\n")

# Huffman编码
huffman_codes, huffman_lengths = huffman_build(p)
huffman_bits = prefix_encode(source, symbols, huffman_codes)
huffman_decoded = prefix_decode(huffman_bits, symbols, huffman_codes)

huffman_average_length = sum(p .* huffman_lengths)
huffman_efficiency = H / huffman_average_length
huffman_compression = 1.0 - length(huffman_bits) / fixed_total_bits

println("【Huffman编码码表】")
print_codebook(symbols, p, huffman_codes)
println("编码码流：", join(huffman_bits))
println("译码序列：", join(huffman_decoded, " "))
println("译码正确：", source == huffman_decoded)
println("总码长：", length(huffman_bits), " bit")
@printf("平均码长：%.6f bit/符号\n", huffman_average_length)
@printf("编码效率：%.4f%%\n", huffman_efficiency * 100)
@printf("相对定长编码压缩率：%.4f%%\n\n", huffman_compression * 100)

# Shannon编码
shannon_codes, shannon_lengths = shannon_build(p)
shannon_bits = prefix_encode(source, symbols, shannon_codes)
shannon_decoded = prefix_decode(shannon_bits, symbols, shannon_codes)

shannon_average_length = sum(p .* shannon_lengths)
shannon_efficiency = H / shannon_average_length
shannon_compression = 1.0 - length(shannon_bits) / fixed_total_bits

println("【Shannon编码码表】")
print_codebook(symbols, p, shannon_codes)
println("编码码流：", join(shannon_bits))
println("译码序列：", join(shannon_decoded, " "))
println("译码正确：", source == shannon_decoded)
println("总码长：", length(shannon_bits), " bit")
@printf("平均码长：%.6f bit/符号\n", shannon_average_length)
@printf("编码效率：%.4f%%\n", shannon_efficiency * 100)
@printf("相对定长编码压缩率：%.4f%%\n\n", shannon_compression * 100)

# 算术编码
arithmetic_bits, final_low, final_high, arithmetic_value =
    arithmetic_encode(source, symbols, p)

arithmetic_decoded =
    arithmetic_decode(arithmetic_bits, length(source), symbols, p)

arithmetic_average_length = length(arithmetic_bits) / length(source)
arithmetic_efficiency = H / arithmetic_average_length
arithmetic_compression = 1.0 - length(arithmetic_bits) / fixed_total_bits

println("【算术编码】")
@printf("最终区间：[%.15f, %.15f)\n", final_low, final_high)
@printf("区间内编码值：%.15f\n", arithmetic_value)
println("编码码流：0.", join(arithmetic_bits))
println("译码序列：", join(arithmetic_decoded, " "))
println("译码正确：", source == arithmetic_decoded)
println("总码长：", length(arithmetic_bits), " bit")
@printf("平均码长：%.6f bit/符号\n", arithmetic_average_length)
@printf("编码效率：%.4f%%\n", arithmetic_efficiency * 100)
@printf("相对定长编码压缩率：%.4f%%\n\n", arithmetic_compression * 100)

println("==================== 性能对比 ====================")
println("编码方法\t总码长\t平均码长\t编码效率\t压缩率")
@printf("Huffman\t%d\t%.4f\t\t%.4f%%\t%.4f%%\n",
        length(huffman_bits), huffman_average_length,
        huffman_efficiency * 100, huffman_compression * 100)
@printf("Shannon\t%d\t%.4f\t\t%.4f%%\t%.4f%%\n",
        length(shannon_bits), shannon_average_length,
        shannon_efficiency * 100, shannon_compression * 100)
@printf("算术编码\t%d\t%.4f\t\t%.4f%%\t%.4f%%\n",
        length(arithmetic_bits), arithmetic_average_length,
        arithmetic_efficiency * 100, arithmetic_compression * 100)
println("==================================================")
