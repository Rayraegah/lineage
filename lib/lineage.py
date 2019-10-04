"""
Lineage – A Genetic Algorithm for Breeding Creatures

@version 1.0.0
@author Rayraegah <rayraegah@gmail.com>

There are two part to the algorithm:
(1) parent gene swapping and
(2) mutation.

(1) Parent gene swapping
    This process is used to get recessive genes into the dominant slot.

    Here's an example:
        There's a creature with traits (D, R1, R2, R3).

    In the above example,
        D = Dominant trait
        R1 = Recessive trait 1
        R2 = Recessive trait 2
        R3 = Recessive trait 3

    When parent gene swapping there's a 25% chance it will swap R3 and R2,
    then R2 with R1, then R1 with D. This happens consequently which means
    the algorithm could end up swapping R3 down to R2, then to R1, then to
    D all in the same iteration. This is extremely unlikely which a chace of
    1.5%, which means R3 genes will less likely to become dominant in one
    step.

(2) Mutation
    This process can only happen after the genes are swapped. The gene swap
    randomizes the genes in different places. Mutation only happens on the
    dominant gene. Once mutated the genes are swapped so that the higher
    valued gene is G2. Some preconditions can be run before a mutation
    like gene2 - gene1 == 1 and is gene1 even or not.

    This algorithm needs both a specific combination of genes and the value
    to be even. The kai values have 1 added to to their binary. All of the
    math should be done in binary.

    After preconditions the mutation happens with a 25% probability if the
    value is less than 23 and half of that otherwise. The mutation is just as
    predictable e.g. mutation = (gene1 / 2) + 16.
"""
import sha3


# get 5-bit chunks of matron and sire
def masker(arg, start, numbytes):
    mask = 2 ** numbytes - 1
    mask = mask << start

    out = arg & mask
    out = out >> start

    return out


# load arguments into bytes arrays in big-Endian order
def get_bytes(genes):
    byte_array = []

    for n in range(32):
        byte_array.append(genes // ((1 << 8) ** n) & 0xff)

    byte_array.reverse()
    return bytes(byte_array)


def get_genemasks(genes):
    genemasks = []

    for n in range(0x30):
        genemasks.append(masker(genes, 5 * n, 5))

    return genemasks


def unmask(masked):
    outs = 0

    for n in range(0x30):
        outs |= masked[n] << 5 * n

    return outs


def breed(matron_genes, sire_genes, matron_cdeb, blockhash):
    # concatenate bytes arrays
    alls = get_bytes(blockhash) + \
        get_bytes(matron_genes) + \
        get_bytes(sire_genes) + \
        get_bytes(matron_cdeb)

    # get hash of bytes arrays. This is your source of "randomness"
    hash = sha3.keccak_256(alls)
    hash = int.from_bytes(hash.digest(), byteorder='big')

    matron_genesmasks = get_genemasks(matron_genes)
    sire_genesmasks = get_genemasks(sire_genes)

    matron_genesmasks_copy = matron_genesmasks.copy()
    sire_genesmasks_copy = sire_genesmasks.copy()

    # note in worst case hashindex wont reach 256 so no need for modulo
    hashindex = 0

    # swap dominant/recessive genes according to masked_hash
    for bigcounter in range(0x0c):
        for smallcounter in range(3, 0, -1):
            count = 4 * bigcounter + smallcounter

            masked_hash = masker(hash, hashindex, 2)
            hashindex += 2
            if masked_hash == 0:
                tmp = matron_genesmasks_copy[count - 1]
                matron_genesmasks_copy[count -
                                       1] = matron_genesmasks_copy[count]
                matron_genesmasks_copy[count] = tmp

            masked_hash = masker(hash, hashindex, 2)
            hashindex += 2
            if masked_hash == 0:
                tmp = sire_genesmasks_copy[count - 1]
                sire_genesmasks_copy[count - 1] = sire_genesmasks_copy[count]
                sire_genesmasks_copy[count] = tmp

    # combine genes from swapped parent genes, introducing mutations
    outmasks = []
    for cnt in range(0x30):
        rando_byte = 0

        # mutate only on dominant genes
        if cnt % 4 == 0:
            tmp1 = matron_genesmasks_copy[cnt] & 1
            tmp2 = sire_genesmasks_copy[cnt] & 1

            if tmp1 != tmp2:
                masked_hash = masker(hash, hashindex, 3)
                hashindex += 3

                mask1 = matron_genesmasks_copy[cnt]
                mask2 = sire_genesmasks_copy[cnt]

                # mutate only if the two parent dominant genes differ by 1...
                if abs(mask2 - mask1) == 1:
                    min_mask = min(mask1, mask2)
                    # and the smaller of the two is even...
                    if min_mask % 2 == 0:
                        if min_mask < 0x17:
                            trial = masked_hash > 1
                        else:
                            trial = masked_hash > 0
                        if not trial:
                            # mutation is the smaller of
                            # the two parent dominant genes,
                            # divided by two, plus 16
                            rando_byte = (min_mask >> 1) + 0x10

            if rando_byte > 0:
                print(cnt)
                outmasks.append(rando_byte)
                continue

        masked_hash = masker(hash, hashindex, 1)
        hashindex += 1

        if masked_hash == 0:
            outmasks.append(matron_genesmasks_copy[cnt])
        else:
            outmasks.append(sire_genesmasks_copy[cnt])

    return unmask(outmasks)
