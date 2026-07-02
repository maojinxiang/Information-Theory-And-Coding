# source_coding_functions.jl
# Huffman、Shannon与算术编码的基础函数
# 适用于基于Julia语法的MWORKS数值计算环境

function check_probability(p)
    if isempty(p)
        error("概率向量不能为空。")
    end
    if any(p .<= 0.0)
        error("所有概率都必须大于0。")
    end
    if abs(sum(p) - 1.0) > 1.0e-10
        error("概率之和必须等于1。")
    end
end

function huffman_build(p)
    check_probability(p)

    n = length(p)
    codes = [Int[] for _ in 1:n]

    if n == 1
        codes[1] = [0]
        return codes, [1]
    end

    probs = copy(p)
    groups = [[i] for i in 1:n]

    while length(probs) > 1
        # 概率升序排列；概率相同时按节点内最小符号序号排序，
        # 使程序每次运行得到稳定一致的码表。
        order = sortperm(1:length(probs),
                         by = i -> (probs[i], minimum(groups[i])))

        probs = probs[order]
        groups = groups[order]

        group0 = groups[1]
        group1 = groups[2]

        for index in group0
            codes[index] = vcat([0], codes[index])
        end

        for index in group1
            codes[index] = vcat([1], codes[index])
        end

        merged_prob = probs[1] + probs[2]
        merged_group = vcat(group0, group1)

        remaining_probs = length(probs) > 2 ? probs[3:end] : Float64[]
        remaining_groups = length(groups) > 2 ? groups[3:end] : Vector{Int}[]

        probs = vcat(remaining_probs, [merged_prob])
        groups = vcat(remaining_groups, [merged_group])
    end

    lengths = length.(codes)
    return codes, lengths
end

function shannon_build(p)
    check_probability(p)

    n = length(p)
    order = sortperm(1:n, by = i -> (-p[i], i))
    p_sorted = p[order]

    codes = [Int[] for _ in 1:n]
    lengths = zeros(Int, n)
    cumulative_probability = 0.0

    for i in 1:n
        code_length = max(1, ceil(Int, -log2(p_sorted[i])))
        x = cumulative_probability
        code = zeros(Int, code_length)

        for k in 1:code_length
            x = 2.0 * x
            bit = floor(Int, x + 1.0e-12)
            bit = clamp(bit, 0, 1)
            code[k] = bit
            x = x - bit
        end

        original_index = order[i]
        codes[original_index] = code
        lengths[original_index] = code_length
        cumulative_probability += p_sorted[i]
    end

    return codes, lengths
end

function prefix_encode(source, symbols, codes)
    bits = Int[]

    for value in source
        index = findfirst(==(value), symbols)

        if isnothing(index)
            error("待编码序列中含有码表未定义的符号。")
        end

        append!(bits, codes[index])
    end

    return bits
end

function prefix_decode(bits, symbols, codes)
    source = eltype(symbols)[]
    position = 1

    while position <= length(bits)
        matched = false

        for i in eachindex(symbols)
            code_length = length(codes[i])
            last_position = position + code_length - 1

            if last_position <= length(bits)
                current_bits = bits[position:last_position]

                if current_bits == codes[i]
                    push!(source, symbols[i])
                    position = last_position + 1
                    matched = true
                    break
                end
            end
        end

        if !matched
            error("码流与当前码表不匹配，无法继续译码。")
        end
    end

    return source
end

function integer_to_bits(integer_value::Int64, bit_number::Int)
    bits = zeros(Int, bit_number)
    value = integer_value

    for i in bit_number:-1:1
        bits[i] = Int(mod(value, 2))
        value = div(value, 2)
    end

    return bits
end

function arithmetic_encode(source, symbols, p)
    check_probability(p)

    cumulative = vcat([0.0], cumsum(p))
    low = 0.0
    high = 1.0

    for value in source
        index = findfirst(==(value), symbols)

        if isnothing(index)
            error("待编码序列中含有概率模型未定义的符号。")
        end

        old_low = low
        range_value = high - low

        low = old_low + range_value * cumulative[index]
        high = old_low + range_value * cumulative[index + 1]

        if high <= low
            error("区间发生精度塌缩，请缩短信源序列。")
        end
    end

    interval_width = high - low
    bit_number = max(1, ceil(Int, -log2(interval_width)) + 1)

    while bit_number <= 52
        scale = Int64(1) << bit_number
        integer_value = ceil(Int64, low * Float64(scale))
        code_value = integer_value / Float64(scale)

        if code_value >= low && code_value < high
            bits = integer_to_bits(integer_value, bit_number)
            return bits, low, high, code_value
        end

        bit_number += 1
    end

    error("未能在双精度有效范围内生成终止码。")
end

function arithmetic_decode(bits, source_length, symbols, p)
    check_probability(p)

    code_value = 0.0
    for i in eachindex(bits)
        code_value += bits[i] * 2.0^(-i)
    end

    cumulative = vcat([0.0], cumsum(p))
    low = 0.0
    high = 1.0
    source = Vector{eltype(symbols)}(undef, source_length)

    for k in 1:source_length
        range_value = high - low

        if range_value <= 0.0
            error("译码区间无效。")
        end

        normalized_value = (code_value - low) / range_value
        normalized_value = clamp(normalized_value, 0.0, 1.0 - 1.0e-12)

        index = length(symbols)
        for i in eachindex(symbols)
            if normalized_value < cumulative[i + 1]
                index = i
                break
            end
        end

        source[k] = symbols[index]

        old_low = low
        low = old_low + range_value * cumulative[index]
        high = old_low + range_value * cumulative[index + 1]
    end

    return source
end

function print_codebook(symbols, p, codes)
    println("符号\t概率\t码字\t码长")
    for i in eachindex(symbols)
        println("$(symbols[i])\t$(round(p[i], digits=4))\t$(join(codes[i]))\t$(length(codes[i]))")
    end
end
