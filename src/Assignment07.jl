module Assignment07

export normalizeDNA,
       composition,
       gc_content,
       complement,
       reverse_complement,
       parse_fasta

# # uncomment the following line if you intend to use BioSequences types
# using BioSequences

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end

function composition(sequence)
    sequence = normalizeDNA(sequence)

    counts = Dict(b=>0 for b in collect("ATGCN"))
    for base in sequence
        counts[base] += 1
    end

    return counts
end

function gc_content(seq)
    comp = composition(seq)
    gc= (comp['G'] + comp['C'] )/ length(seq)
    return gc
end

function complement(base::Char)
    complements = Dict('A'=>'T',
                    'T'=>'A',
                    'G'=>'C',
                    'C'=>'G',
                    'N'=>'N')
    c= complements[uppercase(base)]
    return c
end


function complement(base::AbstractString)
    return join(map(complement, collect(base)))
end

reverse_complement(seq) = reverse(complement(seq))

function parse_fasta(path)
    headers = String[]
    seqs = String[]
    curseq = String[]
    for line in eachline(path)
        if startswith(line, ">")
            push!(headers, line[2:end])
            if length(curseq) > 0
                push!(seqs, join(curseq))
                curseq = String[]
            end
        else
            line = normalizeDNA(line)
            push!(curseq, line)
        end
    end
    length(curseq) > 0 && push!(seqs, join(curseq))
    return headers, seqs
end



end # module Assignment07